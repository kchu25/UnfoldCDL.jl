function get_gt_test(motif::single_block_motif, data)
    mode_parts = [(1,1)];
    lens = [motif.len];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N_test
        if data.raw_data_test[n].mode > 0
            start_ = data.raw_data_test[n].motif_where[1];
            mode_ = data.raw_data_test[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt_test(motif::mixture_k_block_motif, data)
    mode_parts = [(i,1) for i = 1:motif.num_modes];
    lens = [i[end]-i[1]+1 for i in motif.modes];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N_test
        if data.raw_data_test[n].mode > 0
            start_ = data.raw_data_test[n].motif_where[1];
            mode_ = data.raw_data_test[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt_test(motif::gapped_k_block_motif, data)
    mode_parts = [(1,k) for k = 1:motif.K]
    lens = [i.len for i in motif.P_motifs.P]
    lens_rv = reverse(lens)
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N_test
        if data.raw_data_test[n].mode > 0
            start_ = data.raw_data_test[n].motif_where[1];
            end_ = data.raw_data_test[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data_test[n].str[motif_where];
            upper_on = false; 
            pos = [];
            cur_motif = 1;
            counter = 0;
            len_here = data.raw_data_test[n].complement ? lens_rv : lens;
            
            for (ind,s) in enumerate(str)        
                if isuppercase(s) 
                    if !upper_on && counter == 0   
                        push!(pos, start_+ind-1); upper_on = true;                    
                    end      
                    counter+=1;          
                    if counter == len_here[cur_motif]
                        cur_motif += 1;
                        upper_on = false;
                        counter = 0;
                    end                
                else # if it's smaller case
                    upper_on = false;
                end
            end
            if data.raw_data_test[n].complement
                for (ind,p) in enumerate(reverse(pos))
                    covering[mode_parts[ind]][n] = p;
                end
            else
                for (ind,p) in enumerate(pos)
                    covering[mode_parts[ind]][n] = p;
                end
            end                
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt_test(motif::mixture_gapped_k_block_motif, data)
    # use 0 to disregard the mode index.. 
    # TODO: maybe make it more readable later
    mode_parts = [(0,k) for k = 1:motif.motif.K];
    lens = [i.len for i in motif.motif.P_motifs.P];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N_test
       
        if data.raw_data_test[n].mode > 0
            # println(n)
            start_ = data.raw_data_test[n].motif_where[1];
            end_ = data.raw_data_test[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data_test[n].str[motif_where];

            mode_ = data.raw_data_test[n].mode;
            mode_UnitRange = motif.modes[mode_];

            upper_on = false; 
            pos = [];
            cur_motif = 1;
            counter = 0;
            motif_mode = data.raw_data_test[n].complement ? reverse(motif.modes[mode_]) : motif.modes[mode_]
            len_here = data.raw_data_test[n].complement ? lens[motif_mode] : lens[motif_mode]

            for (ind,s) in enumerate(str)                
                if isuppercase(s) 
                    if !upper_on && counter == 0 
                        push!(pos, start_+ind-1); upper_on = true;                    
                    end      
                    counter+=1;       
                    if counter == len_here[cur_motif]         
                        cur_motif += 1;
                        upper_on = false;
                        counter = 0;
                    end                
                else # if it's smaller case
                    upper_on = false;
                end
            end        
            
            # note that length(pos) == length(mode_UnitRange) is true
            @assert length(pos) == length(mode_UnitRange) "$n: $pos $mode_UnitRange"
            if data.raw_data_test[n].complement
                for (ind,p) in enumerate(reverse(pos))
                    covering[(0,mode_UnitRange[ind])][n] = p;
                end
            else
                for (ind,p) in enumerate(pos)
                    covering[(0,mode_UnitRange[ind])][n] = p;
                end
            end        
        end
    end
    ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_overlaps_test(gt::ground_truth, ms, data)
    L = data.L;
    num_fs = length(ms.pfms);
    gt_overlaps = zeros((gt.num_mp, num_fs)); # true positive 
    bg_overlaps = zeros(num_fs);
    activated_f_area = 0;
    activated_f_individial = zeros(num_fs);
    gt_area = 0;
    gt_area_each = zeros(gt.num_mp);
    bg_area = 0;

    # convention: if no key (at nth string) for pfm f, then return -1
    @inbounds for n = 1:data.N_test
        if data.raw_data_test[n].mode > 0
            gt_covering_at_n = covering_at(gt, n);
            sum_g_len = 0.0;
            region_gs = Vector{UnitRange}();
            activated_fs = [Vector{UnitRange}() for _ = 1:ms.num_motifs];

            for g = 1:gt.num_mp 
                g_s = gt_covering_at_n[g]; # if absent, this is -1
                g_e = g_s != -1 ? g_s+gt.lens[g]-1 : -1;
                push!(region_gs, g_s:g_e);
                g_area = g_e-g_s+1;
                gt_area_each[g] += g_s != -1 ? g_area : 0;
                sum_g_len += g_s != -1 ? g_area : 0;
            end
            for f = 1:num_fs
                if haskey(ms.positions[f], n) 
                    for ind = 1:length(ms.positions[f][n])
                        f_s = ms.positions[f][n][ind]
                        f_e = f_s + ms.lens[f]-1;
                        push!(activated_fs[f],  f_s:f_e);
                    end
                end
            end

            region_bg = get_bg_UnitRanges(1:L, region_gs);
            overlap_gt_vs_pfm!(gt_overlaps, region_gs, activated_fs);

            overlap_bg_vs_pfm!(bg_overlaps, region_bg, activated_fs);       
            
            sum_active_vec = sum_range_vec(activated_fs)
            activated_f_individial = activated_f_individial .+ sum_active_vec;
            activated_f_area += sum(sum_active_vec);
            gt_area += sum_g_len; 
            bg_area += L-sum_g_len;
        end
    end
    return gt_area, 
           gt_area_each, 
           bg_area, 
           activated_f_area, 
           activated_f_individial, 
           gt_overlaps, 
           bg_overlaps
end

function get_performance_test(data, ms)
    # note that this modifies the positions and the scores
    overlapping_scan!(ms, data; test=true)
    non_overlap_scan!(ms, data; test=true)
    pvalues = get_fisher_p_values(ms, data; test=true)

    if length(ms.pfms) == 0
        return 0.0,0.0,0.0
    end
    gt = get_gt_test(data.motif, data)
    gt_area, gt_area_each, _, 
    activated_f_area, activated_f_individial, 
        gt_overlaps, bg_overlaps = get_overlaps_test(gt, ms, data);

    sum_overlap = sum(gt_overlaps);
    # gt_cover_ratio = sum(gt_overlaps) ./ gt_area;
    # gt_n_cover = 1 - gt_cover_ratio;

    perf_coeff  = round(sum_overlap/(activated_f_area + gt_area - sum_overlap), digits=2)
    sensitivity = round(sum_overlap/activated_f_area, digits=2)
    specifity   = round(sum_overlap/gt_area, digits=2)
    # fp_ratio = sum(bg_overlaps)/activated_f_area;
    
    return perf_coeff, sensitivity, specifity, pvalues
end
