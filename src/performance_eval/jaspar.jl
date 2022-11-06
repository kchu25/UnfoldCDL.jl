pfm2pwm(pfm) = log2.(pfm ./ .25);
pfm2pwm(pfm, bg) = log2.(pfm ./ bg);

function allr(p, q, p_count, q_count, bg)
    allr_score = Float64[];
    for i = 1:size(p,2)
        view_p_col = Base.view(p, :, i);
        view_q_col = Base.view(q, :, i);
        nb_p = p_count .* view_p_col;
        nb_q = q_count .* view_q_col;
        a1=sum(nb_p .* pfm2pwm(view_q_col, bg)); 
        a2=sum(nb_q .* pfm2pwm(view_p_col, bg));
        push!(allr_score, (a1+a2)/(sum(nb_p)+sum(nb_q)))
    end
    return sum(allr_score)
end

allr_thresh(pfm_len::Integer, fraction_of_columns, score_each_col) = floor(fraction_of_columns*pfm_len)*score_each_col;

function convolve_allr(pfm_c2, pfm,
                       counts_pfm_c2::Integer,
                       counts_pfm::Integer,                        
                       len_pfm_c2::Integer,
                       len_pfm::Integer,
                       thresh_allr::Real,
                       min_col::Integer,
                       bg                   
                       )
    #= len_pfm_c2 will always be smaller since we've select the ones
        with minimal length
    =#;
    allrs = Float64[];
    s1e1s = UnitRange{Int}[];
    s2e2s = UnitRange{Int}[];
    l_dec_1 = Int[]; l_dec_2 = Int[]; 
    r_inc_1 = Int[]; r_inc_2 = Int[];

    for i = 1:(len_pfm_c2+len_pfm-1)
        s1 = max(1, len_pfm_c2-i+1); e1 = min(len_pfm_c2, len_pfm_c2-(i-len_pfm));
        s2 = max(1, i-len_pfm_c2+1); e2 = min(i, len_pfm);
        overlap_count = e1-s1+1;
        push!(s1e1s, s1:e1); push!(s2e2s, s2:e2);
        #=  
            Note that:
            1) no need to calculate if the number of columns of the 
            pfm is less than min_col as specified
            2) no need to calculate the placements for which
            maximal value of the score is below the threshold
        =#
        if overlap_count ≥ min_col && 2*overlap_count ≥ thresh_allr 
            push!(allrs, allr(Base.view(pfm_c2,:,s1:e1), Base.view(pfm,:,s2:e2), 
                            counts_pfm_c2, counts_pfm, bg));
        else
            push!(allrs, -Inf);
        end        
        push!(l_dec_1, max(s2-1,0)); push!(l_dec_2, max(s1-1,0));
        push!(r_inc_1, max(0,len_pfm-i)); push!(r_inc_2, max(i-e2,0));
    end
    argmax_ind = argmax(allrs);
    return allrs[argmax_ind], 
           l_dec_1[argmax_ind], 
           r_inc_1[argmax_ind], 
           l_dec_2[argmax_ind], 
           r_inc_2[argmax_ind]
end

function get_pfm_from_transfac(transfac_path::String; float_type=Float32)
    f = open(transfac_path, "r"); r = readlines(f); close(f);
    rows = Vector{String}();
    start = false
    for i in r
        if i[1:2] == "01" 
            start = true
        elseif i[1:2] == "XX"
            start = false
        end
        start && push!(rows, i);    
    end
    parse_counts = [parse.(float_type, i[2:end]) for i in split.(rows, "\t")];
    count_mat = reduce(hcat, parse_counts); count_mat = count_mat .+ 0.01;
    return count_mat ./ sum(count_mat, dims=1)
end

function get_complement_vec(top10, ms, sort_perm, pfm_jaspar)
    complement_or_not = Bool[];
    for pfm in ms.pfms[sort_perm][1:top10]
        allr_score,_,_,_,_ = Cdlunroll.convolve_allr(pfm_jaspar, pfm, 1000, 1000, size(pfm_jaspar,2), size(pfm,2), 0f0, 3, [.25,.25,.25,.25]);
        allr_score_c,_,_,_,_ = Cdlunroll.convolve_allr(pfm_jaspar, reverse(pfm), 1000, 1000, size(pfm_jaspar,2), size(pfm,2), 0f0, 3, [.25,.25,.25,.25]);
        if allr_score > allr_score_c
            push!(complement_or_not, false)
        else
            push!(complement_or_not, true)
        end            
    end
    return complement_or_not
end

get_result_tuple(matrix_name, 
                 ref_logo_where,
                 top5_pvalues,
                 top5_logo_link,
                 details_link
                 ) = (name=matrix_name, 
                     jaspar_link="https://jaspar.genereg.net/matrix/"*matrix_name, 
                     jaspar_logo=ref_logo_where, 
                     top5_pvalues=top5_pvalues, 
                     top5_logo_link=top5_logo_link, 
                     details_link=details_link
                     )

function get_result_tuple_for_all_jaspar_render(ms, 
                                                pval_str_vec, 
                                                sort_perm,
                                                source_folder_logo_transfac,
                                                ref_save_where,
                                                ref_logo_where,
                                                matrix_name)
    # assume pvalue is a vector and is sorted it in ascending order
    top5 = min(length(pval_str_vec), 5);
    pfm_jaspar = Cdlunroll.get_pfm_from_transfac(source_folder_logo_transfac; 
                        float_type=eltype(ms.pfms[1])); # compare with jaspar logo
    complement_or_not = get_complement_vec(top5, ms, sort_perm, pfm_jaspar)
    top5_pvalues = pval_str_vec[1:top5];
    top5_logo_link = [complement_or_not[i] ? 
                        ref_save_where*"/logos/d$(i)_c.png" : 
                        ref_save_where*"/logos/d$i.png" 
                        for i = 1:top5];
    details_link = ref_save_where*"/summary.html";
    return get_result_tuple(matrix_name, 
                 ref_logo_where,
                 top5_pvalues,
                 top5_logo_link,
                 details_link
                 ) 
end

function save_jaspar_and_get_result_tuple(ms, data,
                                   target_folder, 
                                   source_folder_logo_transfac,
                                   ref_save_where,
                                   ref_logo_where,
                                   matrix_name,
                                   pval_cutoff
                                   )
    isnothing(ms) && return get_result_tuple(matrix_name, ref_logo_where, String[], String[], "");
    
    pvalue_vec, sort_perm = save_result(ms, data, target_folder, pval_cutoff);
    return get_result_tuple_for_all_jaspar_render(ms, 
                                                  pvalue_vec, 
                                                  sort_perm,
                                                  source_folder_logo_transfac,
                                                  ref_save_where,
                                                  ref_logo_where,
                                                  matrix_name)
end
