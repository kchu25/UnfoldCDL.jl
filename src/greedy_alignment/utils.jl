promote_i(x...) = Int.(x);

mutable struct motifs{T <: Integer, S <: Real}
    pfms::Vector{Matrix{S}}
    pwms::Union{Nothing,Vector{Matrix{S}}}
    thresh::Union{Nothing, Vector{S}}
    lens::Vector{T}
    num_motifs::T
    positions::Union{Nothing, Vector{Dict{T, Vector{T}}}}
    scores::Union{Nothing, Vector{Dict{T, Vector{S}}}}
    use_comp::Union{Nothing, Vector{Dict{T, Vector{Bool}}}}
    positions_bg::Union{Nothing, Vector{Dict{T, Vector{T}}}}
    scores_bg::Union{Nothing, Vector{Dict{T, Vector{S}}}}
    use_comp_bg::Union{Nothing, Vector{Dict{T, Vector{Bool}}}}
end

mat_complement(mat::Matrix{T}) where T <: Real = reverse(mat);
get_thresh(pwms, pval_thresh) = float_type.([(pvalue2score(pwm, pval_thresh)) for pwm in pwms]);

pval_thresh_dict_test = Dict{Int, float_type}(6=> 0.00025,
                                              7=> 0.00025,
                                              8=> 0.00025,
                                              9=> 0.00025,
                                              10=>0.00025,
                                              11=>0.00025,
                                              12=>0.00025,
                                              13=>0.00025,
                                              14=>0.00025,
                                              15=>0.00025,
                                              )


pval_thresh_dict = Dict{Int, float_type}(6=> 0.001275,
                                         7=> 0.000925,
                                         8=> 0.0009,
                                         9=> 0.00035,
                                         10=>0.00025,
                                         11=>0.0000895,
                                         12=>0.0000675,
                                         13=>0.00005,
                                         14=>0.0000325,
                                         15=>0.00002,
                                         )

function get_thresh(pwms; test=false)
    lens_info = [size(pwms[i],2) ≥ 15 ? 15 : (size(pwms[i], 2) ≤ 6 ? 6 : size(pwms[i],2)) for i = 1:lastindex(pwms)]
    pval_threshs = ifelse(test, 
                        [pval_thresh_dict_test[l] for l in lens_info], 
                        [pval_thresh_dict[l] for l in lens_info])
    float_type.([(pvalue2score(pwm, pval_threshs[ind])) for (ind,pwm) in enumerate(pwms)]);
end


function mat_compare2!(cmats, cmats_c, ccsum, ccsum_c, i, j, alpha_merge, marked_inds, merged_pair)

        col_k1_leq_col_k2_, offset_, pval_ = 
            compare_k1k2(cmats[i], cmats[j], ccsum[i], ccsum[j])
        col_k1_leq_col_k2_c, offset_c, pval_c = 
            compare_k1k2(cmats[i], cmats_c[j], ccsum[i], ccsum_c[j])
        use_comp = pval_c > pval_;
  
        pval = use_comp ? pval_c : pval_;
        offset = use_comp ? offset_c : offset_;
        col_k1_leq_col_k2 = use_comp ? col_k1_leq_col_k2_c : col_k1_leq_col_k2_;

        if pval > alpha_merge
            marked_inds[i] = true; marked_inds[j] = true;
            push!(merged_pair, (i,j,use_comp, offset, col_k1_leq_col_k2))
        end

end

function merge_count_mats2(count_mats, alpha_merge, diff_tol)
    num_count_mats = length(count_mats)
    lens = size.(count_mats,2)
    marked_inds = fill(false, length(lens))
    merged_pair = Vector{Tuple{Int,Int,Bool,Int,Bool}}()
    # 1st cmat, 2nd cmat, use_comp, offset, 1st cmat smaller or equal to 2nd cmat's length
    count_mats_c = mat_complement.(count_mats)
    count_mat_col_sums = sum.(count_mats, dims=1);
    count_mat_col_sums_c = sum.(count_mats_c, dims=1);
    for i = 1:num_count_mats, j = i+1:num_count_mats
        if  !(marked_inds[i] && marked_inds[j])
            if abs(lens[i]-lens[j]) ≤ diff_tol                
                mat_compare2!(count_mats, count_mats_c,
                        count_mat_col_sums, count_mat_col_sums_c,
                        i, j, alpha_merge, marked_inds, merged_pair)
            end
        end
    end
    new_count_mats = Vector{eltype(count_mats)}(undef, length(merged_pair))
    for (ind, (i, j, use_comp, offset, c1lc2)) in enumerate(merged_pair)
        if c1lc2
            s=offset+1; e=s+size(count_mats[i],2)-1;
            new_count_mats[ind]=count_mats[i]+ifelse(use_comp, (@view count_mats_c[j][:, s:e]), (@view count_mats[j][:, s:e]))
        else
            s=offset+1; e=s+size(count_mats[j],2)-1;
            display((round.(count_mats[i][:,s:e], digits=0)))
            if use_comp
                display(round.(count_mats_c[j], digits=0))
                println("false i and comp")
            else
                display(round.(count_mats[j], digits=0))
                println("false i and not comp")
            end
            println("----")
            new_count_mats[ind]=(@view count_mats[i][:,s:e])+ifelse(use_comp, count_mats_c[j], count_mats[j]);
        end
    end
    vcat(count_mats[setdiff(1:end, findall(marked_inds .== true))], new_count_mats)
end

function mat_compare!(cmats, cmats_c, i, j, alpha_merge, marked_inds, merged_pair)
    if !(marked_inds[i] && marked_inds[j])
        mat2_c = false; is_similar = false;
        if pvalue_g_test(cmats[i], cmats[j]) > alpha_merge
            is_similar = true                    
        elseif pvalue_g_test(cmats[i], cmats_c[j]) > alpha_merge
            is_similar = true
            mat2_c = true
        end
        if is_similar 
            # println((i,j,mat2_c))
            push!(merged_pair, (i,j,mat2_c))
            marked_inds[i] = true; marked_inds[j] = true;
        end
    end
end

# check the pfms of the same length
function merge_count_mats(count_mats, alpha_merge, diff_tol)
    count_mats_c = mat_complement.(count_mats)
    lens = size.(count_mats,2)
    marked_inds = fill(false, length(lens))
    merged_pair = Vector{Tuple{Int,Int,Bool}}()
    uniq_lens = unique(lens)
    for l in uniq_lens
        inds = findall(lens .== l)
        for i = 1:lastindex(inds)
            for j = i+1:lastindex(inds)
                mat_compare!(count_mats, count_mats_c, 
                    inds[i], inds[j], alpha_merge, 
                    marked_inds, merged_pair)
            end
        end
    end
    
    new_count_mats = Vector{eltype(count_mats)}(undef, length(merged_pair))
    for (ind,(ind_1, ind_2, ind_2_comp)) in enumerate(merged_pair)
        new_count_mats[ind] = 
            count_mats[ind_1] + ifelse(ind_2_comp, count_mats_c[ind_2], count_mats[ind_2])
    end
    return vcat(count_mats[setdiff(1:end, findall(marked_inds .== true))], new_count_mats)
end

function merging_count_mats(count_mats, alpha_merge, diff_tol)
    num_count_mats_last = length(count_mats)
    while true
        count_mats = merge_count_mats2(count_mats, alpha_merge, diff_tol)
        num_count_mats = length(count_mats)
        num_count_mats == num_count_mats_last && break
        num_count_mats_last = num_count_mats
    end
    return count_mats
end

function countmats2motifs(count_mats; pval_thresh = 0.00027)
    pfms = (countmat2pfm.(count_mats));
    pwms = freq2pwm.(pfms, float_type)
    score_thresh = get_thresh(pwms)
    lens = size.(pwms,2)
    num_pfms = length(lens)

    return motifs(pfms,
                    pwms,
                    score_thresh,
                    lens,
                    num_pfms,
                    nothing, 
                    nothing,
                    nothing,
                    nothing,
                    nothing,
                    nothing
                );
end

function submat_comlement(data_matrix, start_, end_, k, len)
    piece = @view data_matrix[start_:end_,1,k];
    reverse(reshape(piece,(4,len)))
end

function get_start_end(positions, len4, i, k, ind_k)
    start_ = (positions[i][k][ind_k]-1)*4+1;
    end_ = start_+len4-1;
    return start_, end_
end

function msa_add!(i, len4, positions, use_complement, 
        lens, msa, data_matrix; ps=0.01, return_count_mat=false)
    for k in keys(positions[i])
        for ind_k in 1:length(positions[i][k])
            start_, end_ = get_start_end(positions, len4, i, k, ind_k)
            if use_complement[i][k][ind_k]
                msa .+= 
                    reshape(submat_comlement(data_matrix,start_,end_,k, lens[i]),(4, lens[i]));
            else
                msa .+= reshape((@view data_matrix[start_:end_,1,k]), (4, lens[i]));
            end
        end
    end
    return return_count_mat ? msa .+ ps : countmat2pfm(msa .+ ps)
end

function posdicts2countmats(ms::motifs, keep, data_matrix::Array{S,3}) where {S<:Real}
    count_mats = Vector{Matrix{S}}();
    @inbounds for (i, keep_this) in enumerate(keep)
        if keep_this
            msa = zeros(S, (4, ms.lens[i]));
            count_mat = msa_add!(i, 4*ms.lens[i], ms.positions, 
                ms.use_comp, ms.lens, msa, data_matrix; return_count_mat=true)
            push!(count_mats, count_mat);
        end
    end
    return count_mats
end

function posdict2pos_keeponly(ms, keep::BitVector)
    keep_inds = findall(keep .> 0);
    len_keep = sum(keep);
    pos = [Vector{Tuple{UnitRange{Int64}, Int64, Bool}}() for _ = 1:len_keep]
    q = 1;
    for i in keep_inds
        for k in keys(ms.positions[i])
            for (ind,w) in enumerate(ms.positions[i][k])
                push!(pos[q], (w:w+ms.lens[i]-1, k, ms.use_comp[i][k][ind]))
            end
        end
        q+=1
    end
    return pos
end

