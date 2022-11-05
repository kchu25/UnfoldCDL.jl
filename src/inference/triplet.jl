
const vec_tup_t = NamedTuple{(:at, :f, :len, :f1_prev, :f2_prev, :f1_prev_loc, :f2_prev_loc), 
                        NTuple{7, Int}}

const pos_t = Tuple{UnitRange{Int}, Int};
const topo_t = NamedTuple{(:f1, :f2, :f3, :d12, :d13), NTuple{5, Int}};


const ndict_c_t = Dict{Tuple{Int,Int,Int,Int}, Bool}; # ndict's nested (child) dict type
const ndict_t = Dict{Int, ndict_c_t};

get_code_non_zero_components(Z) = Z[1:4:end,:,:];

function percentile_thresh_Z_nz(Z_nz; percentile= 0.0055)
    Z_nz_sorted = sort((@view Z_nz[:]))
    ind = Int(round(length(Z_nz_sorted) * (1.0 - percentile)))
    return Z_nz_sorted[ind]
end

sort_fcn(x) = x[1];

push_arr!(q, arrs, len) = push!(arrs[q[3]], (at=q[1],
                                                f=q[2], 
                                                len=len, 
                                                f1_prev=0,
                                                f2_prev=0,
                                                f1_prev_loc=0,
                                                f2_prev_loc=0));

function record_condition(i, j; triplet_condition=false)
    if triplet_condition
        cond2_11 = (i.f1_prev != 0 && i.f2_prev != 0) && (j.f1_prev == 0 && j.f2_prev == 0);
        cond2_12 = (i.f1_prev == 0 && i.f2_prev == 0) && (j.f1_prev != 0 && j.f2_prev != 0);
        cond2_1 = cond2_11 || cond2_12;
        cond2_2 = true
        ((i.f1_prev == j.f) && (i.f1_prev_loc == j.at)) && (cond2_2=false);
        ((i.f2_prev == j.f) && (i.f2_prev_loc == j.at)) && (cond2_2=false);
        ((j.f1_prev == i.f) && (j.f1_prev_loc == i.at)) && (cond2_2=false);
        ((j.f2_prev == i.f) && (j.f2_prev_loc == i.at)) && (cond2_2=false);
        return cond2_1 && cond2_2;
    else
        return true
    end
end

function record_pairs_n!(pair_dist_record, arr_n_i, arr_n_j, max_start_apart; triplet_condition=false)
    if record_condition(arr_n_i, arr_n_j; triplet_condition=triplet_condition)
        k1, k2, pos_diff = arr_n_i.f, arr_n_j.f, arr_n_j.at - arr_n_i.at
        @inbounds if 1 ≤ pos_diff ≤ max_start_apart
            pair_dist_record[k1,k2,pos_diff] += 1;
        end
    end
end

function record_pairs!(pair_dist_record, arrs_n, max_start_apart; triplet_condition=triplet_condition)
    this_len = length(arrs_n);
    @inbounds for i = 1:this_len, j = i+1:this_len
        record_pairs_n!(pair_dist_record, arrs_n[i], arrs_n[j], max_start_apart; triplet_condition=triplet_condition)
    end
end

#=
arrs: array that contain arrays of filters and their activate positions
base_flen: the filter length in the learning
max_flen: maximum filter length of all the filters
M: number of filters (more if triplet_condition = true)
seq_len: sequence length in the dataset (should be all the same)
=#
function get_code_pair_enrichment(arrs, max_flen, M, N; triplet_condition=false, max_start_apart_incr = 3)
    max_start_apart = 2*max_flen + max_start_apart_incr;
    pair_dist_record = zeros(UInt16, (M,M,max_start_apart)); 
    # UInt16: largest value is 2^16 -1
    # https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/
    @info "filling the code pair enrichment..."
    @inbounds for n = 1:N
        n % 500 == 0 && (@info "Done $(round(100*(n/N), digits=1)) %")
        record_pairs!(pair_dist_record, arrs[n], max_start_apart; triplet_condition=triplet_condition)
    end
    return pair_dist_record
end

function get_filter_activations_sorted_array(code, flen, thresh)
    found = findall(code .> thresh);
    arrs = [Vector{vec_tup_t}() for _ = 1:size(code,3)];
    for q in found push_arr!(q, arrs, flen); end
    sort!.(arrs , by=sort_fcn);
    return arrs
end

function pair_cond_1(code, l, f1, f2, n , Lc, l_p_d, thresh)
    code[l,f1,n] > thresh && l_p_d ≤ Lc && code[l_p_d,f2,n] > thresh;
end

function get_arr_pair_i!(i, ind, Lc, M, N, arrs2, code, thresh, max_flen, flen; nextto=false)
    f1, f2, d = nextto ? (i[1], i[2], flen) : (i[1], i[2], i[3]);
    for n = 1:N, l = 1:Lc                          
        l_p_d = l+d;
        if pair_cond_1(code, l, f1, f2, n, Lc, l_p_d, thresh) 
            flen_here = d+flen; 
            max_flen[1] = flen_here > max_flen[1] ? flen_here : max_flen[1];
            push!(arrs2[n], (at=l, f=M+ind, len=flen_here, 
                             f1_prev=f1, f2_prev=f2, f1_prev_loc=l, f2_prev_loc=l_p_d))
        end
    end
end

function get_pair_filter_activations_sorted_array(_lap, code, thresh, count_threshold, flen; nextto=false)
    Lc, M, N = size(code);    
    arrs2 = [Vector{vec_tup_t}() for _ = 1:N];
    max_flen = [0];
    found2 = findall(_lap .> count_threshold);
    num_pairs = length(found2);
    for (ind,i) in enumerate(found2)
        get_arr_pair_i!(i, ind, Lc, M, N, arrs2, code, thresh, max_flen, flen; nextto=nextto)
    end
    return arrs2, max_flen[1], num_pairs;
end

function push_arr!(arr1_n, arr2_n)
    for i in arr2_n
        push!(arr1_n, i)
    end
end

function get_arrs_prev(arrs, fil_num)
    arrs_ = [Vector{vec_tup_t}() for _ = 1:data.N];
    for (n,i) in enumerate(arrs)
        for j in i
            j.f == fil_num && (push!(arrs_[n], j));
        end
    end
    return arrs_
end

function find_f1f2_prev(arrs2, fil_num)
    break_cond = false;
    f1_prev = nothing; f2_prev = nothing; d = nothing;
    for i in arrs2
        for j in i
            if j.f == fil_num
                break_cond=true;
                f1_prev = j.f1_prev;
                f2_prev = j.f2_prev;
                d = j.f2_prev_loc - j.f1_prev_loc;
            end
            break_cond && break;
        end
        break_cond && break;
    end
    @assert !isnothing(f1_prev) "f1_prev not found"
    @assert d ≥ 1 "d should be larger than or equal to 1"
    return f1_prev, f2_prev, d
end

function get_three(f1,f2,d,M,arrs2)
    # get the triples and the distance between (f1,f2) and (f1,f3)
    # determine whether fp1 fp2 are at the front
    # either f1 > M or f2 > M
    # @assert !(f1 > M && f2 > M) "Only one of filters must belong to the pairs"
    @assert (f1 ≤ M) ⊻ (f2 ≤ M) " $f1 $f2 Only one of filters must belong to the pair/nonpair"
    fp1fp2_front = f1 > M; 
    fp1, fp2, dp12 = fp1fp2_front ? find_f1f2_prev(arrs2, f1) : find_f1f2_prev(arrs2, f2);
    f_1, f_2, f_3 = fp1fp2_front ? (fp1, fp2, f2) : (f1, fp1, fp2);
    d_12, d_13 = fp1fp2_front ? (dp12, d) : (d, d+dp12); # important !
    return f_1, f_2, f_3, d_12, d_13
end

function triple_cond(l2, l3, Lc, l, f_1, f_2, f_3, n, code, thresh)
    l2 ≤ Lc && l3 ≤ Lc && code[l,f_1,n] > thresh && code[l2,f_2,n] > thresh && code[l3,f_3,n] > thresh
end

function chk_triple(ref_dict, f_1, f_2, f_3, d_12, d_13)
    # check if this triple is something we've recorded before
    # if f_1 == 73 && f_2 == 28 && f_3 == 63 && d_12== 2 && d_13 == 3
    #     println(f_1, " ",  f_2, " ", f_3, " ", d_12, " ", d_13)
    #     println("(d12: $d_12, d13: $d_13)")
    #     println(haskey(ref_dict[f_1], (f_2,f_3)))
    #     println(all(ref_dict[f_1][f_2,f_3] .== (d_12, d_13)))
    #     println("ref_dict[f_1][f_2,f_3]:  $(ref_dict[f_1][f_2,f_3])")
    #     println(haskey(ref_dict[f_1], (f_3,f_2)))
    #     println(all(ref_dict[f_1][f_3,f_2] .== (d_13, d_12)))
    #     println("ref_dict[f_1][f_3,f_2]:  $(ref_dict[f_1][f_3,f_2])")
    # end

    @inbounds if haskey(ref_dict, f_1)
        if haskey(ref_dict[f_1], (f_2,f_3,d_12,d_13))
            return false
        end
        if haskey(ref_dict[f_1], (f_3,f_2,d_13,d_12))
            return false
            # if all(ref_dict[f_1][f_3,f_2] .== (d_13, d_12))
            # if ref_dict[f_1][f_3,f_2][1] == d_13 &&  ref_dict[f_1][f_3,f_2][2] == d_12

                # return false
            # end
        end
        return true        
    end
    return true
end

function add_ref!(ref_dict, f_1, f_2, f_3, d_12, d_13)
    # add reference to the reference dict
    !haskey(ref_dict, f_1) && (ref_dict[f_1] = ndict_c_t())        
    ref_dict[f_1][f_2, f_3, d_12, d_13] = true;
end

triplet_cover_len(d_12, d_13, flen) = max(d_12, d_13) + flen;

function get_pos_i(N, Lc, f_1, f_2, f_3, d_12, d_13, code, thresh, flen)
    pos_i = pos_t[];
    for n = 1:N, l = 1:Lc
        l2 = l+d_12;
        l3 = l+d_13;
        if triple_cond(l2, l3, Lc, l, f_1, f_2, f_3, n, code, thresh)                
            push!(pos_i, (l:l+triplet_cover_len(d_12,d_13,flen)-1, n))
        end
    end
    return pos_i
end

function get_triple_pos_i!(i, arrs2, Lc, M, N, flen, code, thresh, pos, ref_dict; nexto=false)
    f1, f2, d = nexto ? (i[1], i[2], flen) : (i[1], i[2], i[3]);
    f_1, f_2, f_3, d_12, d_13 = get_three(f1, f2, d, M, arrs2)    

    if d_12 != d_13 && chk_triple(ref_dict, f_1, f_2, f_3, d_12, d_13)

        push!(pos, get_pos_i(N, Lc, f_1, f_2, f_3, d_12, d_13, code, thresh, flen));

        add_ref!(ref_dict, f_1, f_2, f_3, d_12, d_13); # for ref dict 

        return (f1=f_1, f2=f_2, f3=f_3, d12=d_12, d13=d_13);
    else
        return nothing
    end
end

function get_positions_triples(_lap2, count_threshold, arrs2, code, thresh, flen; nexto=false)
    Lc, M, N = size(code)
    found = findall(_lap2 .> count_threshold); num_found = length(found);
    pos = Vector{Vector{pos_t}}();
    topologies = topo_t[];

    ref_dict = ndict_t();  # reference dictionary for already calculated triplets
    @info "Getting the triples..."
    println("# found: $(length(found))")
    for (ind,i) in enumerate(found)
        ind % 500 == 0 && (@info "Done $(round(100*(ind/num_found), digits=1)) %")
        topology = get_triple_pos_i!(i, arrs2, Lc, M, N, flen, code, thresh, pos, ref_dict; nexto=nexto);
        !isnothing(topology) && push!(topologies, topology);
    end
    pos, topologies
end

function get_enriched_pairs(Z_nz, flen, thresh)
    # Z_nz: non-zero components of the code (inside the constraint)
    arrs = get_filter_activations_sorted_array(Z_nz, flen, thresh);
    # sort!.(arrs , by=sort_fcn);
    pair_dist_record = get_code_pair_enrichment(arrs, flen, size(Z_nz,2), size(Z_nz,3));
    return arrs, pair_dist_record 
end

function get_triples(_lap, 
                     arrs, 
                     code, 
                     thresh, 
                     pair_count_thresh, 
                     flen; 
                     nexto=false)

    arrs_cp = copy(arrs);
    arrs2, max_flen, num_pairs = 
        get_pair_filter_activations_sorted_array(_lap, code, thresh, pair_count_thresh, flen; nextto=nexto);
    push_arr!.(arrs_cp, arrs2);
    sort!.(arrs_cp, by=sort_fcn);
    pair_dist_record = 
        get_code_pair_enrichment(arrs_cp, max_flen, size(code,2)+num_pairs, size(code,3); triplet_condition=true); 
        # extra space needed here, as seen by size(code,2)+num_pairs -- may do some data structure improvement later
    return arrs2, pair_dist_record
end

function get_pos_topo(_lap2, triple_count_threshold, arrs2, code, thresh, flen; nexto=false)
    pos, topologies = 
        get_positions_triples(_lap2, triple_count_threshold, arrs2, code, thresh, flen; nexto=nexto);
    return pos, topologies
end

function get_position_topo_all(Z_nz, flen, thresh; 
                               pair_count_thresh=20, 
                               triple_count_thresh=20,
                               count_thresh=10, 
                               r=100,         # upperbound for number of outliers
                               alpha=0.05     # significance level
                               )
    # Z_nz: the (nonzero components) sparse code returned from the CDL training
    # flen: the filter len
    arrs, pair_dist_record = 
        get_enriched_pairs(Z_nz, flen, thresh);
    arrs2, triple_dist_record = 
        get_triples(pair_dist_record, arrs, Z_nz, thresh, pair_count_thresh, flen);    
    pos, topologies = 
        get_pos_topo(triple_dist_record, triple_count_thresh, arrs2, Z_nz, thresh, flen);
    num_outliers = get_number_of_outliers(triple_dist_record; 
        count_thresh=count_thresh, r=r, alpha=alpha) # generalized esd test
    
    sort_inds = sortperm(pos, by=length, rev=true);
    return pos[sort_inds][1:min(length(pos),num_outliers)], topologies[sort_inds][1:min(length(topologies),num_outliers)]
end

################## for testing synthetic data ########################

function num_overlap(r1::UnitRange, r2::UnitRange)
    @assert r1[1] ≤ r1[end] && r2[1] ≤ r2[end] "range is not valid"
    if r1[1] ≤ r2[1] ≤ r1[end]
        return min(r1[end],r2[end])-r2[1]+1;
    elseif  r2[1] ≤ r1[1] ≤ r2[end]
        return min(r1[end],r2[end])-r1[1]+1;
    else
        return 0;
    end
end

function check_syn_data_olap(pos, topologies, data, flen, num)
    olap2_record = [];
    for (ind, p) in enumerate(pos)
        this_cover_len = triplet_cover_len(topologies[ind].d12, topologies[ind].d13, flen)
        overlap = 0;
        for q in p
            overlap += num_overlap(data.raw_data[q[2]].motif_where, q[1]);
        end
        push!(olap2_record, (overlap / (length(p)*this_cover_len), length(p)))
    end
    # olap2_record
    @info "Mean cover percent: $(mean([i[1] for i in olap2_record]))"

    sort_cover = sortperm([i[1] for i in olap2_record], rev=true);

    # println(sort_cover[1])
    for i = 1:num
        @info "Max cover $i: $(olap2_record[sort_cover[i]])"
    end
    for i = 1:num
        @info "Max by triplet $(topologies[sort_cover[i]])"
    end
    return sort_cover
end

function chk_olap_perc_i(pos, topologies, flen, data, i)
    this_cover_len = triplet_cover_len(topologies[i].d12, topologies[i].d13, flen)
    overlap = 0;
    for q in pos[i]
        overlap += num_overlap(data.raw_data[q[2]].motif_where, q[1]);
    end
    overlap / (length(pos[i])*this_cover_len)
end


function get_counts_and_coverperc(pos, topologies, data, flen)
    counts = Int[];
    perc = Float64[];
    for (ind, p) in enumerate(pos)
        this_cover_len = triplet_cover_len(topologies[ind].d12, topologies[ind].d13, flen)
        overlap = 0;
        for q in p
            overlap += num_overlap(data.raw_data[q[2]].motif_where, q[1]);
        end
        # push!(olap2_record, (overlap / (length(p)*this_cover_len), length(p)))
        push!(counts, length(p))
        push!(perc, overlap / (length(p)*this_cover_len))
    end
    sorted_ind = sortperm(counts);
    return counts[sorted_ind], perc[sorted_ind]
end

######################################################################

