function scan_trainset_cpu!(ms, data, dep::Bool)
    overlapping_scan!(ms, data);
    dep && non_overlap_scan!(ms, data)
    overlapping_scan_bg!(ms, data)
    non_overlapping_scan_bg!(ms, data)
end

function scan_trainset_gpu!(ms, data, dep::Bool)
    scan_w_gpu!(ms, data)
    dep && non_overlap_scan!(ms, data)
    scan_w_gpu!(ms, data; scan_bg=true)
    non_overlapping_scan_bg!(ms, data)
end

function scan_testset!(ms, data)
    overlapping_scan!(ms, data; test=true);
    non_overlap_scan!(ms, data; test=true); 
    overlapping_scan_bg!(ms, data; test=true)
    non_overlapping_scan_bg!(ms, data; test=true)
end

function run(ms, data, alpha_fisher, ic_expand_t, ic_shrink_t, diff_tol, allr_thresh, gpu; dep=false)
    gpu ? scan_trainset_gpu!(ms, data, dep) : scan_trainset_cpu!(ms, data, dep);
    pvals = get_fisher_p_values_(ms, data);
    keep_f = pvals .< alpha_fisher
    keep_e = expansions_ms!(ms, data, keep_f; ic_expand_t=ic_expand_t, ic_shrink_t=ic_shrink_t)
    # get rid of small PWMs (less than 4 columns)
    keep = keep_e .& keep_f .& (ms.lens .> 4)
    count_mats_ = posdicts2countmats(ms, keep, data.data_matrix);
    if dep
        ms_pos = posdict2pos_keeponly(ms, keep_f)
        c1 = merging_count_mats(count_mats_, ms_pos, data, ic_shrink_t, diff_tol; allr_thresh=allr_thresh)
        return countmats2motifs(c1);
    else
        return countmats2motifs(count_mats_);
    end
end

data_matrix_gpu(data) = 
    cu(reshape(data.data_matrix, (size(data.data_matrix,1), data.N)));

get_data_matrix_bg_gpu(data) = 
    cu(reshape(data.data_matrix_bg, (size(data.data_matrix_bg,1), data.N_test)));

function greedy_alignment(new_cmats, data;
                          ic_expand_t=0.7,
                          ic_shrink_t=0.4,
                          alpha_fisher = 1e-5,
                          allr_thresh =0.25,
                          indep_run=2,
                          dep_run=5,
                          diff_tol=4,
                          gpu=false)
    ms = countmats2motifs(new_cmats);
    data.data_matrix_gpu = reshape(data.data_matrix_gpu, (size(data.data_matrix_gpu,1), data.N));
    data.data_matrix_bg_gpu = reshape(data.data_matrix_bg_gpu, (size(data.data_matrix_bg_gpu,1), data.N));

    for _ = 1:indep_run
        ms = run(ms, data, alpha_fisher, ic_expand_t, ic_shrink_t, diff_tol, allr_thresh, gpu)        
    end
    for _ = 1:dep_run
        ms = run(ms, data, alpha_fisher, ic_expand_t, ic_shrink_t, diff_tol, allr_thresh, gpu; dep=true)
    end
    return ms
end