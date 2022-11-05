"""
    Z: the sparse code returned from the training
    data: dataset
    flen: the length of the filters (they are all of the same length)
"""
function obtain_count_matrices(Z,                                  
                                data, 
                                flen;
                                percentage_thresh=0.0055,
                                k=6,
                                esd_alpha=0.05,
                                esd_r=1000,
                                esd_count_thresh=10,
                                ic_expand_t=0.7,
                                ic_shrink_t=0.4,
                                allr_thresh=0.25,
                                diff_tol=4
                                )
    # get the enriched triplets
    Z_nz = get_code_non_zero_components(Z);
    thresh = percentile_thresh_Z_nz(Z_nz; percentile=percentage_thresh)
    pos = obtain_pos(Z_nz, thresh, flen; k=k, count_thresh=esd_count_thresh, r=esd_r, alpha=esd_alpha)
    new_pos = get_expanded_pos(pos, data; 
                               ic_expand_thresh=ic_expand_t, 
                               ic_shrink_thresh=ic_shrink_t
                               );

    sort!(new_pos, by=length, rev=true)
    count_mats = countmat.([get_strs(data, new_pos, k) for k = 1:length(new_pos)]);
    new_cmats = merging_count_mats(count_mats, new_pos, data, ic_shrink_t, diff_tol; allr_thresh=allr_thresh)
    
    return new_cmats
end