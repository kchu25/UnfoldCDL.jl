# greedy search kernel
function greedy_search!(pwms, data_dat_gpu, lens, pos_scores, thresh)
    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # kth pair
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y; # nth sequence
    l = (blockIdx().z - 1) * blockDim().z + threadIdx().z; # lth position

    L, N = size(data_dat_gpu); L_div_4 = CUDA.Int(L/4);
    K, _, _ = size(pwms);
    if k ≤ K && n ≤ N && l ≤ L_div_4-lens[k]+1
        @inbounds for (ind,i) in enumerate(l:l+lens[k]-1)
            for a = 1:4
                pos_scores[k,n,l] += pwms[k,a,ind]*data_dat_gpu[(i-1)*4+a,n];                
            end
        end
        pos_scores[k,n,l] = pos_scores[k,n,l] > thresh[k] ? pos_scores[k,n,l] : 0f0;
    end
    return nothing
end

function modify_w_found!(found, positions, scores, use_comp, pos_scores_arr; rc=false)
    comp = rc ? true : false;
    @inbounds for f in found
        m, n, l = f[1], f[2], f[3];
        if haskey(positions[m], n)
            push!(positions[m][n], l)
            push!(scores[m][n],    pos_scores_arr[m,n,l])
            push!(use_comp[m][n],  comp)
        else
            positions[m][n] = [l];
            scores[m][n]    = [pos_scores_arr[m,n,l]];
            use_comp[m][n]  = [comp];
        end
    end
end

function get_pos_scores_arr(ms, data; scan_bg=false, rc=false, S=Float32)
    data_matrix_gpu = scan_bg ? data.data_matrix_bg_gpu : data.data_matrix_gpu;

    maxlen = maximum(ms.lens);
    pwms = zeros(S, ms.num_motifs, 4, maxlen);
    @inbounds for i = 1:ms.num_motifs pwms[i,:,1:ms.lens[i]] = 
        rc ? reverse(ms.pwms[i]) : ms.pwms[i]; end
    pos_scores = CUDA.zeros(S, ms.num_motifs, data.N, data.L);
    @cuda threads=ker_3d blocks=b_size_3d(pos_scores) greedy_search!(cu(pwms), 
                                                              data_matrix_gpu, 
                                                              cu(ms.lens), 
                                                              pos_scores,
                                                              cu(ms.thresh));
    pos_scores_arr = Array(pos_scores);
    found = findall(pos_scores_arr .> 0);
    return pos_scores_arr, found
end

function gpu_scan(ms, data; scan_bg=false)
    pos_scores_arr, found = get_pos_scores_arr(ms, data; scan_bg=scan_bg, rc=false);
    pos_scores_rc_arr, found_rc = get_pos_scores_arr(ms, data; scan_bg=scan_bg, rc=true);
    
    # motifs prep
    positions, scores, use_comp = motifs_prep(ms);
    modify_w_found!(found, positions, scores, use_comp, pos_scores_arr; rc=false)
    modify_w_found!(found_rc, positions, scores, use_comp, pos_scores_rc_arr; rc=true)
    return positions, scores, use_comp
end

function scan_w_gpu!(ms, data; scan_bg=false)
    
    positions, scores, use_comp = gpu_scan(ms, data; scan_bg=scan_bg)

    if scan_bg
        ms.positions_bg = positions;
        ms.scores_bg = scores;
        # ms.use_comp_bg = use_comp; # not needed
    else
        ms.positions = positions;
        ms.scores = scores;
        ms.use_comp = use_comp;
    end
end
