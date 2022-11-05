struct ground_truth{T <: Integer, T2 <: Integer}
    mode_parts::Vector{Tuple{T,T}}              # (mode, k)
    num_mp::T
    covering::Dict{Tuple{T,T}, Dict{T2,T2}}     # key=(mode,k), val=Dict(n=>pos); val[n] may be empty, and that's fine
    lens::Vector{T}
end

function get_gt(motif::single_block_motif, data)
    mode_parts = [(1,1)];
    lens = [motif.len];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            mode_ = data.raw_data[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt(motif::mixture_k_block_motif, data)
    mode_parts = [(i,1) for i = 1:motif.num_modes];
    lens = [i[end]-i[1]+1 for i in motif.modes];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            mode_ = data.raw_data[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt(motif::gapped_k_block_motif, data)
    mode_parts = [(1,k) for k = 1:motif.K]
    lens = [i.len for i in motif.P_motifs.P]
    lens_rv = reverse(lens)
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            end_ = data.raw_data[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data[n].str[motif_where];
            upper_on = false; 
            pos = [];
            cur_motif = 1;
            counter = 0;
            len_here = data.raw_data[n].complement ? lens_rv : lens;
            
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
            if data.raw_data[n].complement
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

function get_gt(motif::mixture_gapped_k_block_motif, data)
    # use 0 to disregard the mode index.. 
    # TODO: maybe make it more readable later
    mode_parts = [(0,k) for k = 1:motif.motif.K];
    lens = [i.len for i in motif.motif.P_motifs.P];
    covering = Dict(mp=>Dict{Int32,Int32}() for mp in mode_parts);
    for n = 1:data.N
       
        if data.raw_data[n].mode > 0
            # println(n)
            start_ = data.raw_data[n].motif_where[1];
            end_ = data.raw_data[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data[n].str[motif_where];

            mode_ = data.raw_data[n].mode;
            mode_UnitRange = motif.modes[mode_];

            upper_on = false; 
            pos = [];
            cur_motif = 1;
            counter = 0;
            motif_mode = data.raw_data[n].complement ? reverse(motif.modes[mode_]) : motif.modes[mode_]
            len_here = data.raw_data[n].complement ? lens[motif_mode] : lens[motif_mode]

            # len_here = data.raw_data[n].complement ? lens_rv : lens;

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
            if data.raw_data[n].complement
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


#= 
return an ordered (according to mode_parts) array of 
covering starts for each mode_parts in seqeunce n 
=#
function covering_at(g::ground_truth, n::Integer) 
    starts_at = Vector{Integer}();
    # for all the (mode, k) in the ground truth motif
    for (mode,k) in g.mode_parts
        if haskey(g.covering[(mode,k)], n)
            push!(starts_at, g.covering[(mode,k)][n]);
        else
            push!(starts_at, -1);
        end
    end
    return starts_at
end

overlap(r1::UnitRange, r2::UnitRange) = r1[1] ≤ r2[1] ≤ r1[end] || r1[1] ≤ r2[end] ≤ r1[end] || r2[1] ≤ r1[1] ≤ r2[end] || r2[1] ≤ r1[end] ≤ r2[end];
not_found(g_start::Integer, f_start::Integer) = g_start == -1 || f_start == -1;
appeared(r::UnitRange) = r[1] != -1 && r[end] != -1;

# sum_range_vec(rs::Vector{UnitRange}) = 
#     [r[1] != -1 ? r[end]-r[1]+1 : 0 for r in rs];

sum_range_each(r) = r[end]-r[1]+1
# sum_range_vec(rs::Vector{Vector{UnitRange}}) = map(x->sum(sum_range_each.(x)), rs)

function sum_range_vec(rs::Vector{Vector{UnitRange}})
    reduced_rs = reduce(vcat, rs)
    length(reduced_rs) == 0 && return 0
    return sum(sum_range_each.(reduced_rs))
end

sum_range(rs::Vector{UnitRange}) = sum(sum_range_vec(rs));

#= 
Return the unit ranges in bg that does not intersect with ground truths
e.g. r_bg = 1:100
     r_gts = [2:6, 8:15, 55:60]
     return [1:1, 7:7, 16:54, 61:100]     

    test:
    a = 1:100;
    r1 = 2:6;
    r2 = 8:15;
    r3 = 55:60;
    rs = [r1,r2,r3];
    get_bg_UnitRanges(a,rs)
=#
function get_bg_UnitRanges(r_bg::UnitRange, r_gts::Vector{UnitRange}) 
    set_r_gts = Set{Int}();
    for r in r_gts union!(set_r_gts, Set(r)) end
    bg_wo_r_gts = sort(collect(setdiff(Set(r_bg), set_r_gts)));

    last=0; bg = Vector{Vector{Int}}(); push!(bg,Vector{Int}()); counter = 1;
    for i in bg_wo_r_gts        
        if last == i-1 || last == 0
            push!(bg[counter], i)
        else
            push!(bg, Vector{Int}())
            counter += 1;
            push!(bg[counter], i)
        end
        last = i;
    end
    return [i[1]:i[end] for i in bg]
end

function overlap_gt_vs_pfm!(gt_overlap_matrix::Matrix{Q},
                                 r_gts,
                                 r_fs
                                 ) where Q <: Real
    for (g_ind, rg) in enumerate(r_gts)
        for (f_ind, rf) in enumerate(r_fs)
            for this_range in rf
                gt_overlap_matrix[g_ind, f_ind] += num_overlap(rg, this_range);                        
            end
        end
    end
end

#=
increment the overlap count between the background 
and each of the pfms in a sequence 
=#
function overlap_bg_vs_pfm!(overlap_vec::Vector{Q},
                                 r_bgs,       
                                 r_fs                         
                                ) where Q <: Real
    for (_, rb) in enumerate(r_bgs)
        for (f_ind, rf) in enumerate(r_fs)
            for this_range in rf
                overlap_vec[f_ind] += num_overlap(rb,this_range);
            end
        end
    end
end

function get_overlaps(gt::ground_truth, ms, data)
    L = Int(size(data.data_matrix,1)/4);
    num_fs = length(ms.pfms);
    gt_overlaps = zeros((gt.num_mp, num_fs)); # true positive 
    bg_overlaps = zeros(num_fs);
    activated_f_area = 0;
    activated_f_individial = zeros(num_fs);
    gt_area = 0;
    gt_area_each = zeros(gt.num_mp);
    bg_area = 0;

    # convention: if no key (at nth string) for pfm f, then return -1
    @inbounds for n = 1:data.N
        if data.raw_data[n].mode > 0
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
                # f_s = haskey(ms.positions[f], n) ? ms.positions[f][n] : -1;
                # f_e = f_s != -1 ? f_s + ms.lens[f]-1 : -1;
                # push!(activated_fs,  f_s:f_e);
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

function get_performance(data, ms)

    overlapping_scan!(ms, data; test=false)
    non_overlap_scan!(ms, data; test=false)

    if length(ms.pfms) == 0
        return 0.0,0.0,0.0
    end
    gt = get_gt(data.motif, data)
    gt_area, gt_area_each, _, 
    activated_f_area, activated_f_individial, 
        gt_overlaps, bg_overlaps = get_overlaps(gt, ms, data);

    sum_overlap = sum(gt_overlaps);

    perf_coeff  = round(sum_overlap/(activated_f_area + gt_area - sum_overlap), digits=2)
    sensitivity = round(sum_overlap/activated_f_area, digits=2)
    specifity   = round(sum_overlap/gt_area, digits=2)
    
    return perf_coeff, sensitivity, specifity    
end
