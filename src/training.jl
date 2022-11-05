function fill_Z!(Z, i, end_ind, Z_)
    @inbounds Z[:,:,i:end_ind] .= Array(Z_); # TODO: Can we do this faster??
end

function retrieve_sparse_representation(data_matrix, len, projs, hp, ffts, cdl; F=float_type)
    @info "Retreiving the sparse code and the filters..."
    Z = Array{F,3}(undef, (len.C,hp.M,len.Nf));
    ZᵀZ_sum = CUDA.zeros(F, (len.cc_vlen, hp.M, 1));
    ZᵀS_sum = CUDA.zeros(F, (len.cs_vlen, hp.M, 1));

    @inbounds for i = 1:hp.batch_size:len.Nf
        end_ind = i+hp.batch_size-1 ≤ len.Nf ? i+hp.batch_size-1 : len.Nf
        d_gpu = cu(data_matrix[:,:,i:end_ind]);
        Z_ = CSC(d_gpu, cdl, len, hp, projs);
        fill_Z!(Z, i, end_ind, Z_)
        zeros_Zz, zeros_Zs, zeros_Sz, rFzz, irFzz, 
            rFzs, irFzs, rFsz = batch_prep_du(size(d_gpu,3), len, projs, ffts);
        Zr = reverse(Z_, dims=1);
        ZᵀZ, ZᵀS = BregmanPrep(Z_, Zr, d_gpu, rFzz, irFzz, rFzs, rFsz, irFzs, zeros_Zz, zeros_Zs, zeros_Sz)
        ZᵀZ_sum .+= sum(ZᵀZ, dims=3)
        ZᵀS_sum .+= sum(ZᵀS, dims=3)
    end

    ZᵀZ_sum = Array(ZᵀZ_sum); ZᵀS_sum = Array(ZᵀS_sum);
    D_final, mu = square_non_neg_params_du(cdl)
    projs_mapfrange_cpu = Array(projs.mapfrange);
    D_cpu = Array(D_final);
    ZᵀZ_sum = reshape(ZᵀZ_sum, (len.cc_vlen, 1, hp.M));
    DZᵀZ    = convolution(ZᵀZ_sum, D_cpu, pad=len.f_len_inc);
    DZᵀZ_sum = sum(DZᵀZ, dims=2); 
    D_grad   = reshape(
        batched_mul(projs_mapfrange_cpu,  ZᵀS_sum + reshape(DZᵀZ_sum, (len.cs_vlen,hp.M,1))),
            (size(projs.mapfrange,1), 1, hp.M));
    Breg_num = reshape(D_cpu .* exp.(-mu[1] .* D_grad), (4, hp.filter_len, 1, hp.M));
    D_final_cpu = reshape((Breg_num ./ sum(Breg_num,dims=1)), (hp.f_len, 1, hp.M))
    return D_final_cpu, Z
end


function feed_dataloader(data_matrix, hp)
    return Flux.DataLoader(data_matrix, batchsize=hp.batch_size, shuffle=true)
end

function training_prep(M, f, batch_size, K_c, gamma_1, pool_mask_len, pool_stride, data_matrix, learning_rate)
    hp              = Cdlunroll.HyperParam(batch_size, f, f*4, M, K_c, 1, gamma_1, pool_mask_len, pool_stride);
    data_load       = Cdlunroll.feed_dataloader(data_matrix, hp);
    len             = Cdlunroll.length_info(data_load, data_matrix, hp)
    projs           = Cdlunroll.projectors(hp, len);
    ffts            = Cdlunroll.cuda_fft_plans(len, hp);
    cdl             = Cdlunroll.CDL(hp, 1, K_c+1);
    ps              = Flux.params(cdl);
    # opt             = Flux.AdaBelief()
    # opt             = Flux.AdaBelief(0.0002)
    opt             = Flux.AdaBelief(learning_rate)
    return hp, data_load, len, projs, ffts, cdl, ps, opt
end

function train_basic(data_matrix;
                     M             = 65, 
                     num_epochs    = 25, 
                     f             = 7, 
                     K_c           = 3, 
                     gamma_1       = 23.5, 
                     batch_size    = 120,
                     pool_mask_len = 40,
                     pool_stride   = 4,
                     learning_rate = 0.0003
                     )

    N = size(data_matrix,3);   
    batch_size = N < batch_size ? N-1 : batch_size;
    @info "# filters: $M , filter length: $f, epochs: $num_epochs, # coding layers: $K_c, γ: $gamma_1, batch size: $batch_size, learning_rate: $learning_rate"
    N % batch_size == 0 && (batch_size = batch_size + 1)

    hp, data_load, len, projs, ffts, cdl, ps, opt =
        training_prep(M, f, batch_size, K_c, gamma_1, 
            pool_mask_len, pool_stride, data_matrix, learning_rate)

    for iter = 1:num_epochs
        for d in data_load
            d = d |> gpu;
            gs = gradient(ps) do 
                CDLforward(d, cdl, len, hp, projs, ffts)
            end
           
            Flux.Optimise.update!(opt, ps, gs) ## update parameters                
        end

        # just to show the training loss
        if iter % 5 == 0 
            @info "$iter epoch done."
            # @info "μ: $(cdl.mu)"
            # @info "λ: $(cdl.lambda)"
            loss_value = CDLforward(first(data_load) |> gpu, cdl, len, hp, projs, ffts);
            @info "epoch: $iter, batch loss: $loss_value"
            # isnan(loss_value) && (@info "nan loss"; return error("nan gradient"))
        end
    end
    # @info "Obtaining the sparse code..."
    # D, _  = CDLforward(first(data_load) |> gpu, cdl, len, hp, projs, ffts; get_parameters=true);
    # Z = CSC_full(data_load.data |> gpu, cdl, len, hp, projs);
    D, Z = retrieve_sparse_representation(data_matrix, len, projs, hp, ffts, cdl)
    return Z, D
end



