mutable struct length_info
    L::Int          # length of each DNA sequence (times 4) 
    Nf::Int         # number of all training data points 
    N::Int          # number of (batched) data points 
    Ns::Int         # number of (small-batched) data points 
    MN::Int         # (number of filters) times (number of data points)
    f_len_inc::Int  # filter increment length when convolved with filters
    C::Int          # code length
    CN::Int         # code length times number of (batched) data points
    cc_vlen::Int    # convolved length for ZᵀZ 
    ccf_vlen::Int   # convolved length for DZᵀZ
    cs_vlen::Int    # convolved length for ZᵀS  (subject to change for small batch size)
    cf_vlen::Int    # convolved length for ZᵀD
    cff_vlen::Int   # convolved length for DᵀDX    
    fs_vlen::Int    # convolved length for DᵀS  (subject to change for small batch size)
    C_nz::Int       # number of non-zero code components of xₘₙ
    C_pool::Int     # code length after max pooling

    function length_info(data_load_gpu, data_matrix_gpu, hp; wo_label=true)
        L = wo_label ? size(data_load_gpu.data, 1) : size(data_load_gpu.data.data, 1); 
        Nf = size(data_matrix_gpu,3);
        N = hp.batch_size; 
        N_s = wo_label ? size(data_load_gpu.data,3) % N : size(data_load_gpu.data.data,3) % N;
        M = hp.M; 
        f_len_inc = hp.f_len-1;
        C = L-f_len_inc; 
        cc_vlen = 2*C-1;
        ccf_vlen = cc_vlen+f_len_inc;
        cs_vlen = C+L-1;
        cf_vlen = C+f_len_inc;
        cff_vlen = cf_vlen+f_len_inc;
        fs_vlen = L+f_len_inc;
        C_nz = Int(floor(C/4+1));
        C_pool = Int(floor((C - hp.pool_mask_len)/hp.pool_stride)+1);

        new(L, Nf, N, N_s, M*N, f_len_inc, C, C*N, cc_vlen,
            ccf_vlen, cs_vlen, cf_vlen, cff_vlen, fs_vlen, C_nz, C_pool)
    end
end

mutable struct projectors
    zeros_Zz_n::CuArray{float_type, 3}             # (C-1,M,N)
    zeros_Zs_n::CuArray{float_type, 3}             # (L-1,M,N)
    zeros_Sz_n::CuArray{float_type, 3}             # (C-1,1,N)
    zeros_Zz_s::CuArray{float_type, 3}             # (C-1,M,N)
    zeros_Zs_s::CuArray{float_type, 3}             # (L-1,M,N)
    zeros_Sz_s::CuArray{float_type, 3}             # (C-1,1,N)
    zeros_Zz_1::CuArray{float_type, 3}             # (C-1,M,N)
    zeros_Zs_1::CuArray{float_type, 3}             # (L-1,M,N)
    zeros_Sz_1::CuArray{float_type, 3}             # (C-1,1,N)
    mapfrange::CuArray{float_type, 2}
    z_mask_n::CuArray{float_type, 3, CUDA.Mem.DeviceBuffer}
    z_mask_s::CuArray{float_type, 3, CUDA.Mem.DeviceBuffer}
    z_mask_1::CuArray{float_type, 3, CUDA.Mem.DeviceBuffer}
    qI::CuArray{float_type, 2}

    function projectors(hp, len)
        @info "Constructing projectors..."
        # mapfrange   = spzeros(hp.f_len, len.cs_vlen);
        mapfrange   = zeros(hp.f_len, len.cs_vlen);
        mapfrange[1:hp.f_len, len.C:len.C+len.f_len_inc] = Matrix(I, hp.f_len, hp.f_len);
        mapfrange   = cu(mapfrange);
        
        # code mask
        z_mask_col = 1:len.C .∈ [1:4:len.C]; 
        z_mask_n = cu(float_type.(repeat(z_mask_col, outer=(1,hp.M,len.N))));
        z_mask_s = cu(float_type.(repeat(z_mask_col, outer=(1,hp.M,len.Ns))));
        z_mask_1 = cu(float_type.(repeat(z_mask_col, outer=(1,hp.M,1))));

        new(CUDA.zeros(float_type, (len.C-1, hp.M, len.N)),
            CUDA.zeros(float_type, (len.L-1, hp.M, len.N)),
            CUDA.zeros(float_type, (len.C-1, 1,    len.N)),
            CUDA.zeros(float_type, (len.C-1, hp.M, len.Ns)),
            CUDA.zeros(float_type, (len.L-1, hp.M, len.Ns)),
            CUDA.zeros(float_type, (len.C-1, 1,    len.Ns)),
            CUDA.zeros(float_type, (len.C-1, hp.M, 1)),
            CUDA.zeros(float_type, (len.L-1, hp.M, 1)),
            CUDA.zeros(float_type, (len.C-1, 1,    1)),
            mapfrange,
            z_mask_n,
            z_mask_s,
            z_mask_1,
            cu(Matrix(I, hp.M, hp.M))
        )
    end
end

struct cuda_fft_plans
    # for normal batch
    rFzz_n::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzz_n::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}
    rFzs_n::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    rFsz_n::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzs_n::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}

    rFzz_s::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzz_s::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}
    rFzs_s::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    rFsz_s::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzs_s::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}

    rFzz_1::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzz_1::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}
    rFzs_1::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    rFsz_1::CUDA.CUFFT.rCuFFTPlan{float_type, -1, false, 3}
    irFzs_1::AbstractFFTs.ScaledPlan{ComplexF32, CUDA.CUFFT.rCuFFTPlan{ComplexF32, 1, false, 3}, float_type}

    function cuda_fft_plans(len, hp)
        @info "Constructing FFT plans..."
        Zr_padded_z = CuArray{float_type,3}(undef, (len.cc_vlen, hp.M, len.N));
        Zr_padded_s = CuArray{float_type,3}(undef, (len.cs_vlen, hp.M, len.N));
        S_padded  = CuArray{float_type,3}(undef, (len.cs_vlen,1,len.N));

        rFzz_n  = plan_rfft(Zr_padded_z, 1); fzz1 = rFzz_n*Zr_padded_z; 
        irFzz_n = plan_irfft(fzz1, len.cc_vlen, 1);
        rFzs_n  = plan_rfft(Zr_padded_s, 1);
        rFsz_n  = plan_rfft(S_padded, 1); fzs = rFzs_n*Zr_padded_s;
        irFzs_n = plan_irfft(fzs,len.cs_vlen, 1);

        Zr_padded_z = CuArray{float_type,3}(undef, (len.cc_vlen, hp.M, len.Ns));
        Zr_padded_s = CuArray{float_type,3}(undef, (len.cs_vlen, hp.M, len.Ns));
        S_padded  = CuArray{float_type,3}(undef, (len.cs_vlen,1,len.Ns));

        rFzz_s  = plan_rfft(Zr_padded_z, 1); fzz1 = rFzz_s*Zr_padded_z; 
        irFzz_s = plan_irfft(fzz1, len.cc_vlen, 1);
        rFzs_s  = plan_rfft(Zr_padded_s, 1);
        rFsz_s  = plan_rfft(S_padded, 1); fzs = rFzs_s*Zr_padded_s;
        irFzs_s = plan_irfft(fzs,len.cs_vlen, 1);

        Zr_padded_z = CuArray{float_type,3}(undef, (len.cc_vlen, hp.M, 1));
        Zr_padded_s = CuArray{float_type,3}(undef, (len.cs_vlen, hp.M, 1));
        S_padded  = CuArray{float_type,3}(undef, (len.cs_vlen, 1, 1));

        rFzz_1  = plan_rfft(Zr_padded_z, 1); fzz1 = rFzz_1*Zr_padded_z; 
        irFzz_1 = plan_irfft(fzz1, len.cc_vlen, 1);
        rFzs_1  = plan_rfft(Zr_padded_s, 1);
        rFsz_1  = plan_rfft(S_padded, 1); fzs = rFzs_1*Zr_padded_s;
        irFzs_1 = plan_irfft(fzs, len.cs_vlen, 1);

        new(
            rFzz_n, irFzz_n, rFzs_n, rFsz_n, irFzs_n,
            rFzz_s, irFzz_s, rFzs_s, rFsz_s, irFzs_s,
            rFzz_1, irFzz_1, rFzs_1, rFsz_1, irFzs_1
        )
    end
end

