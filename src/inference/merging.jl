get_pos_len(pos) = length(pos[1][1])
seq_num(pos, ind) = pos[ind][2];
front_ind(pos, ind) = pos[ind][1][1];
mode_dict_incr!(mode_dict, diff) =
    haskey(mode_dict, diff) ? mode_dict[diff] += 1 : mode_dict[diff] = 1

function get_max_key(mode_dict)
    # gives the offset value of pos_k1 and pos_k2 that is the mode
    maxkey, maxvalue = first(mode_dict)
    for (key, value) in mode_dict
      value > maxvalue && (maxkey = key; maxvalue = value)
    end
    maxkey
end

ratio_overlap(mode_dict::Dict, offset_mode::Int, pos_ki) = 
    mode_dict[offset_mode]/length(pos_ki)

function find_diff_in_same_seq2!(pos_k1, ind_k1, pos_k2, seq_num_k1, 
                                 mode_dict, used_seq_num)
    @inbounds for ind_k2 in 1:length(pos_k2)
        seq_num_k2 = seq_num(pos_k2, ind_k2);
        if seq_num_k1 == seq_num_k2
            used_seq_num[seq_num_k1] = true;
            mode_dict_incr!(mode_dict, 
                front_ind(pos_k1, ind_k1) - front_ind(pos_k2, ind_k2))
        end
    end
end

function get_mode_dict(pos_k1, pos_k2)
    mode_dict = Dict{Int, Int}();
    used_seq_num = Dict{Int, Bool}();

    for ind_k1 in 1:length(pos_k1)
        seq_num_k1 = seq_num(pos_k1, ind_k1);
        !haskey(used_seq_num, seq_num_k1) && 
        find_diff_in_same_seq2!(pos_k1, ind_k1, pos_k2, seq_num_k1, mode_dict, used_seq_num)
    end
    return mode_dict
end

function pvals_chisq(q; deg_of_freedom=3)
    csq = Chisq(deg_of_freedom);
    pvalues = Array{Float64}(undef, length(q))
    @inbounds for (ind,i) in enumerate(q)
        pvalues[ind] = 1-cdf(csq,i)
    end
    return pvalues
end

function pvalue_g_test(cmat1, cmat2)
    ccsum1, ccsum2 = sum(cmat1, dims=1), sum(cmat2, dims=1);
    pvalue_g_test(cmat1, cmat2, ccsum1, ccsum2)
end

function pvalue_g_test(cmat1, cmat2, ccsum1, ccsum2)
    N_acgt = cmat1 + cmat2
    N      = ccsum1 + ccsum2
    N_jb_1 = (N_acgt .* ccsum1) ./ N
    N_jb_2 = (N_acgt .* ccsum2)  ./ N
    g = sum(2 .* cmat1 .* log.(cmat1 ./ N_jb_1), dims=1) +
        sum(2 .* cmat2 .* log.(cmat2 ./ N_jb_2), dims=1);    
    pval_gtest   = pvals_chisq(g)

    # k = Int(ceil(0.75*size(cmat1,2)));
    # pval_gtest = pval_gtest[sortperm(pval_gtest)[1:k]];
    # println(pval_gtest)
    
    geomean_pval = geomean(pval_gtest);
    mean_pval    = mean(pval_gtest);
    # println("g mean $geomean_pval")
    # println("a mean $mean_pval")
    # ga_mean      = (geomean_pval+mean_pval) ./ 2
    # ga_mean      = (0.2*geomean_pval+0.8*mean_pval)
    return mean_pval
    # return sum(-log.(pval_gtest))
end

function g_test_at_offset(cmat1, cmat2, ccsum1, ccsum2)
    # assume cmat2's number of column is longer than or equal to cmat1
    pvals = Float64[]
    c1_size = size(cmat1,2);
    c2_size = size(cmat2,2);
    @inbounds for i = 1:(c2_size-c1_size+1)
        push!(pvals, 
            pvalue_g_test(cmat1, (@view cmat2[:,i:i+c1_size-1]), 
                          ccsum1, (@view ccsum2[:,i:i+c1_size-1]))
                          )
    end
    return pvals
end

sort_pval_olap(x) = (x[4],x[5])

function compare_k1k2(cmat_k1, cmat_k2, ccsum1, ccsum2)
    pvals = nothing;
    col_k1_leq_col_k2 = size(cmat_k1,2) ≤ size(cmat_k2, 2)
    if col_k1_leq_col_k2
        pvals = g_test_at_offset(cmat_k1, cmat_k2, ccsum1, ccsum2)
    else
        pvals = g_test_at_offset(cmat_k2, cmat_k1, ccsum2, ccsum1)
    end
    max_ind = argmax(pvals);    
    return col_k1_leq_col_k2, max_ind-1, pvals[max_ind]
end

function get_adjusted_pos2(pos_ki, front_Δ, back_Δ, L)
    # adjust the "smaller" positions to be larger
    adjusted = Vector{eltype(pos_ki)}(undef, length(pos_ki))
    todelete = Int[];
    @inbounds for (ind, p) in enumerate(pos_ki)
        start_ = p[1][1]-front_Δ; 
        end_   = p[1][end]+back_Δ;
        if start_ ≥ 1 && end_ ≤ L
            adjusted[ind] = (start_:end_, p[2])
        else
            push!(todelete, ind)
        end
    end
    adjusted[setdiff(1:end, todelete)]
end

function get_adjusted_pos3(pos_ki, front_Δ, back_Δ)
    # adjust the "larger" positions to be smaller
    adjusted = Vector{eltype(pos_ki)}(undef, length(pos_ki))
    @inbounds for (ind, p) in enumerate(pos_ki)
        start_ = p[1][1]+front_Δ; 
        end_   = p[1][end]-back_Δ;
        adjusted[ind] = (start_:end_, p[2])        
    end
    adjusted
end

function merge_position!(offset, pos_ki_small, pos_kj_large, L)
    front_Δ = offset
    back_Δ = get_pos_len(pos_kj_large)-(get_pos_len(pos_ki_small)+offset)
    unique(vcat(pos_kj_large, 
                get_adjusted_pos2(pos_ki_small, front_Δ, back_Δ, L)))
    # unique(vcat(get_adjusted_pos3(pos_kj_large, front_Δ, back_Δ), pos_ki_small))
end

function get_olap_ratio(pos_k1, pos_k2)::Float64
    mode_dict = get_mode_dict(pos_k1, pos_k2)
    isempty(mode_dict) && return 0.0
    offset_mode = get_max_key(mode_dict)
    ratio_overlap(mode_dict, offset_mode, pos_k2)
end

function merge_k1k2!(pos, offset, col_k1_leq_col_k2, k1, k2, merged, merged_pos, L)
    if col_k1_leq_col_k2
        push!(merged_pos, merge_position!(offset, pos[k1], pos[k2], L))
    else
        push!(merged_pos, merge_position!(offset, pos[k2], pos[k1], L))
    end
    merged[k1] = true; 
    merged[k2] = true;
end

function merge_k1_positions2!(k1, num_pos, data, pos, merged, merged_pos, 
                              cmats, cmat_colsums, diff_tol, alpha_merge)
    k1_len = get_pos_len(pos[k1])

    @inbounds for k2 in k1+1:num_pos        
        k2_len = get_pos_len(pos[k2])
        if abs(k1_len-k2_len) ≤ diff_tol
            col_k1_leq_col_k2, offset, pval = 
                compare_k1k2(cmats[k1], cmats[k2], cmat_colsums[k1], cmat_colsums[k2])
            pval > alpha_merge && (merge_k1k2!(pos, offset, col_k1_leq_col_k2, k1, k2, merged, merged_pos, data.L); break)     
        end
    end
end

function merge_position2(pos, data; diff_tol=4, alpha_merge=0.1)
    num_pos = length(pos);
    merged = fill(false, length(pos));
    count_mats = countmat.([get_strs(data, pos, k) for k = 1:num_pos]);
    count_mat_col_sums = sum.(count_mats, dims=1);
    merged_pos = Vector{eltype(pos)}();
    
    @inbounds for k1 in 1:num_pos
        !merged[k1] && merge_k1_positions2!(k1, num_pos, data,
            pos, merged, merged_pos, count_mats, count_mat_col_sums, 
            diff_tol, alpha_merge);
    end
    new_pos = vcat(pos[.!merged], merged_pos)
    sorted_ind = sortperm(new_pos, by=length, rev=true)
    return new_pos[sorted_ind]
end

function allr(countmat₁, countmat₂, pfm₁, pfm₂)
    pwm₁ = log2.(pfm₁ ./ 0.25);
    pwm₂ = log2.(pfm₂ ./ 0.25);
    allr_vec = sum(countmat₂ .* pwm₁ .+ countmat₁ .* pwm₂, dims=1) ./ 
            sum(countmat₁ + countmat₂, dims=1)
    return sum(allr_vec)/length(allr_vec)
end

function allr_at_offset(cmat1, cmat2, pfm1, pfm2)
    # assume cmat2's number of column is longer than or equal to cmat1
    allrs = Float64[]
    c1_size = size(cmat1,2);
    c2_size = size(cmat2,2);
    @inbounds for i = 1:(c2_size-c1_size+1)
        push!(allrs, 
            allr(cmat1, (@view cmat2[:,i:i+c1_size-1]), 
                    pfm1, (@view pfm2[:,i:i+c1_size-1]))
                          )
    end
    return allrs
end

mutable struct mat
    count_mat::Matrix{float_type}
    count_mat_c::Matrix{float_type}
    pfm::Matrix{float_type}
    pfm_c::Matrix{float_type}
    len::Int
    processed::Bool
end

struct merged_info
    ind::Int
    lq_this::Bool # less than or equal to this count matrix
    comp::Bool
    allr::float_type
    len::Int
    offset_i::Union{Int, Nothing}
    offset_j::Union{Int, Nothing}
    right_inc::Union{Int, Nothing}
end

get_offset_i(mi::merged_info) = isnothing(mi.offset_i) ? 0 : mi.offset_i;
get_max_offset_i(vec_info::Vector{merged_info}) = maximum(get_offset_i.(vec_info))
get_right_inc(mi::merged_info) = isnothing(mi.right_inc) ? 0 : mi.right_inc;
get_max_right_inc(vec_info::Vector{merged_info}) = maximum(get_right_inc.(vec_info))
get_j_ind(mi::merged_info) = mi.ind;
get_js_ind(vec_info::Vector{merged_info}) = get_j_ind.(vec_info);

function info(t::merged_info)
    @info "ind: $(t.ind)"
    @info "i_lq_this: $(t.lq_this)"
    @info "comp: $(t.comp)"
    @info "allr: $(t.allr)"
    @info "len: $(t.len)"
    @info "offset_i: $(t.offset_i)"
    @info "offset_j: $(t.offset_j)"
    @info "right_inc $(t.right_inc)"
end

diff_small(t1::mat, t2::mat, diff_tol) = abs(t1.len-t2.len) < diff_tol;
get_count_mat(t::mat) = t.count_mat;

function merging_count_mats(count_mats, new_pos, data, ic_shrink_t, diff_tol; allr_thresh=0.25)
    count_mats_c = Cdlunroll.mat_complement.(count_mats);
    pfms = Cdlunroll.countmat2pfm.(count_mats);
    pfms_c = reverse.(pfms);

    @inbounds my_cmats = [mat(count_mats[i], 
                  count_mats_c[i],
                  pfms[i],
                  pfms_c[i],
                  size(count_mats[i],2), false) for i = 1:lastindex(count_mats)];

    new_cmats = Vector{Matrix{float_type}}();

    for i_1 = 1:lastindex(my_cmats)
        if !my_cmats[i_1].processed       
            my_cmats[i_1].processed = true;
            merged_i_1_info = Vector{merged_info}();        

            for i_2 = 1:lastindex(my_cmats)            
                if !my_cmats[i_2].processed && 
                    diff_small(my_cmats[i_1], my_cmats[i_2], diff_tol)
                    i_lq_j, j_comp, allr, offset_i, offset_j, right_inc = 
                            compare_two_triplets_allr(my_cmats[i_1], my_cmats[i_2]);
                    if allr > allr_thresh
                        my_cmats[i_2].processed = true
                        push!(merged_i_1_info, 
                            merged_info(i_2, i_lq_j, j_comp, allr, 
                            my_cmats[i_2].len, offset_i, offset_j, right_inc));
                    end
                end
            end

            # println(i_1)
            if !isempty(merged_i_1_info)
                this_count_mat = get_result_merged_mat(merged_i_1_info, 
                        new_pos, data, ic_shrink_t, my_cmats[i_1].len, i_1);
                size(this_count_mat,2) > 0 && push!(new_cmats, this_count_mat)
            else
                size(my_cmats[i_1].count_mat,2) > 0 && push!(new_cmats, my_cmats[i_1].count_mat)
            end
        end
    end
    return new_cmats
end

function get_result_merged_mat(merged_i_1_info, 
                               new_pos, 
                               data, 
                               ic_shrink_t, 
                               len_i,
                               i_1; 
                               merge_up=true
                               )
    inds, front_Δ, back_Δ = get_increments(merged_i_1_info, len_i, i_1; merge_up=merge_up)
    modified_vec_positions = all_vec_pos_incr(new_pos, inds, front_Δ, back_Δ);
    merged_count_mat = get_merged_count_mat(data, modified_vec_positions, merged_i_1_info);
    trimmed_count_mat = get_trimmed_count_mat(merged_count_mat , ic_shrink_t);
    return trimmed_count_mat
end

function get_merged_count_mat(data, modified_vec_positions, merged_i_1_info)
    those_count_mats = Cdlunroll.countmat.(
                [Cdlunroll.get_strs_ignore_out_of_range(data, position)
                    for position in modified_vec_positions]);
    merged_count_mat = those_count_mats[1];
    for (ind,cm) in enumerate(@view those_count_mats[2:end])
        merged_count_mat += merged_i_1_info[ind].comp ? reverse(cm) : cm;
    end
    return merged_count_mat
end

function get_j_increments_merge_up(merged_i_1_info::Vector{merged_info}, 
                          len_i, 
                          i_front_Δ, 
                          i_back_Δ
    )
    front_Δ = zeros(Int, length(merged_i_1_info))
    back_Δ  = zeros(Int, length(merged_i_1_info))
    for (ind, mi) in enumerate(merged_i_1_info)
        if !mi.lq_this
            @assert !isnothing(mi.offset_j) "offset_j nothing"
            Δ = len_i - (mi.offset_j+mi.len)
            front_Δ[ind] = -(i_front_Δ + (mi.comp ? Δ : mi.offset_j));
            back_Δ[ind]  = i_back_Δ + (mi.comp ? mi.offset_j : Δ);
        else
            Δ = mi.len - (mi.offset_i+len_i);         

            front_Δ[ind] = -(i_front_Δ - (mi.comp ? Δ : mi.offset_i));
            back_Δ[ind]  = (i_back_Δ - (mi.comp ? mi.offset_i : Δ))
            # println(ind, " ",  mi.len, " ", (mi.offset_i+len_i), " ", len_i)
        end
    end
    return front_Δ, back_Δ
end

function get_j_increments_no_merge_up(merged_i_1_info::Vector{merged_info},   
                            len_i
                        )
    front_Δ = zeros(Int, length(merged_i_1_info))
    back_Δ  = zeros(Int, length(merged_i_1_info))
    # figure out the signs here so that only addition operation is needed later on 
    for (ind, mi) in enumerate(merged_i_1_info)
        if !mi.lq_this
            @assert !isnothing(mi.offset_j) "offset_j nothing"
            Δ = len_i - (mi.offset_j+mi.len)
            front_Δ[ind] = -(mi.comp ? Δ : mi.offset_j);
            back_Δ[ind]  = mi.comp ? mi.offset_j : Δ;
        else
            Δ = mi.len - (mi.offset_i+len_i)
            front_Δ[ind] = mi.comp ? Δ : mi.offset_i;
            back_Δ[ind]  = -(mi.comp ? mi.offset_i : Δ);
        end
    end
    return front_Δ, back_Δ
end

function get_increments(merged_i_1_info::Vector{merged_info}, len_i, i_1; merge_up=true)
    i_front_Δ = merge_up ? get_max_offset_i(merged_i_1_info) : 0;
    i_back_Δ  = merge_up ? get_max_right_inc(merged_i_1_info)  : 0;

    js_front_Δ, js_back_Δ = merge_up ? get_j_increments_merge_up(merged_i_1_info, 
                                   len_i, 
                                   i_front_Δ, 
                                   i_back_Δ
                                   ) : get_j_increments_no_merge_up(merged_i_1_info, len_i);
    js_inds = get_js_ind(merged_i_1_info)
    return pushfirst!(js_inds, i_1),
           pushfirst!(js_front_Δ, -i_front_Δ), 
           pushfirst!(js_back_Δ, i_back_Δ)
end

pos_incr(p::Tuple{UnitRange{Int}, Int, Bool}, f_Δ::Int, b_Δ::Int) = 
    ((p[1][1]+f_Δ):(p[1][end]+b_Δ), p[2], p[3])

vec_pos_incr(ps::Vector{Tuple{UnitRange{Int}, Int, Bool}}, f_Δ::Int, b_Δ::Int) = 
    pos_incr.(ps, f_Δ, b_Δ)

all_vec_pos_incr(
    ps_vec::Vector{Vector{Tuple{UnitRange{Int}, Int, Bool}}}, inds, front_Δ, back_Δ) = 
    vec_pos_incr.((@view ps_vec[inds]), front_Δ, back_Δ)

function ic_col(count_mat) # assume bg = 0.25
    pfm = Cdlunroll.countmat2pfm(count_mat);
    2 .+ sum(pfm .* log2.(pfm), dims=1)
end

get_trim_coor(m, ic_shrink_t) = ic_col(m) .< ic_shrink_t;

function get_trimmed_count_mat(count_mat, ic_shrink_t)
    c = get_trim_coor(count_mat, ic_shrink_t)
    left_start = 1; right_stop=length(c);
    while left_start ≤ length(c) && c[left_start] == 1  left_start+=1 end
    while right_stop ≥ 1 && c[right_stop] == 1 right_stop-=1 end
    return count_mat[:, left_start:right_stop]
end

function compare_two_triplets_allr(cmat_i, cmat_j)
    i_lq_j = cmat_i.len ≤ cmat_j.len 
    if i_lq_j
        allrs_fls_ij = Cdlunroll.allr_at_offset(cmat_i.count_mat, 
                                        cmat_j.count_mat, 
                                        cmat_i.pfm, 
                                        cmat_j.pfm)
        allrs_fls_ijc = Cdlunroll.allr_at_offset(cmat_i.count_mat, 
                                          cmat_j.count_mat_c, 
                                          cmat_i.pfm, 
                                          cmat_j.pfm_c)
        max_allr_ij_ind  = argmax(allrs_fls_ij); 
        max_allr_ijc_ind = argmax(allrs_fls_ijc);

        allr_ij  = allrs_fls_ij[max_allr_ij_ind]
        allr_ijc = allrs_fls_ijc[max_allr_ijc_ind]

        j_comp      = allr_ijc > allr_ij;
        allr_i_lq_j = j_comp ? allr_ijc : allr_ij;
        offset      = j_comp ? max_allr_ijc_ind-1 : max_allr_ij_ind-1;
        
        right_inc = cmat_j.len-(offset+cmat_i.len)

        return i_lq_j, j_comp, allr_i_lq_j, offset, nothing, right_inc
    else
        allrs_fls_ji = Cdlunroll.allr_at_offset(cmat_j.count_mat, 
                                        cmat_i.count_mat, 
                                        cmat_j.pfm, 
                                        cmat_i.pfm)
        allrs_fls_jci = Cdlunroll.allr_at_offset(cmat_j.count_mat_c, 
                                          cmat_i.count_mat, 
                                          cmat_j.pfm_c, 
                                          cmat_i.pfm)
        max_allr_ji_ind  = argmax(allrs_fls_ji); 
        max_allr_jci_ind = argmax(allrs_fls_jci);
        
        allr_ji  = allrs_fls_ji[max_allr_ji_ind]
        allr_jci = allrs_fls_jci[max_allr_jci_ind]

        j_comp     = allr_jci > allr_ji;
        allr_j_l_i = j_comp ? allr_jci : allr_ji
        offset     = j_comp ? max_allr_jci_ind-1 : max_allr_ji_ind-1;

        return i_lq_j, j_comp, allr_j_l_i, nothing, offset, nothing
    end
end