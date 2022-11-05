sort_cartesian(x) = (x[3], x[1]);

function get_n_range(cartesian_array, N)
    # given a cartesian_array, find the range of indices that indicate each sequence
    range_ind_arr=zeros(Int, (N,2)); n=0; num_rows=size(cartesian_array,1);
    for i = 1:num_rows
        if cartesian_array[i][3] != n
            n+=1; range_ind_arr[n,1] = i;
            n > 1 && (range_ind_arr[n-1,2] = i-1);
        end
    end
    range_ind_arr[end,2]=num_rows
    return range_ind_arr
end

function get_indicator_arrays(Z_nz, thresh)
    found = findall(Z_nz .> thresh);
    sort!(found, by=sort_cartesian);
    range_ind_arr = get_n_range(found, size(Z_nz,3));
    return found, range_ind_arr
end

function binary_lower(value_low_bdd, c_arr, n)
    start_ = 1
    end_ = n
    ans = -1
    while(start_ ≤ end_)
        mid = (start_+end_) ÷ 2
        if c_arr[mid][1] ≥ value_low_bdd
            end_ = mid-1
            ans = mid
        else
            start_ = mid+1
        end
    end
    return ans
end

function get_valid_Triplet_ranges(c_arr, k)
    ranges = Tuple{Int,Int}[]; n = length(c_arr)
    for i = 3:n
        # smallest index that has value ≥ arr[i]-k
        si = binary_lower(c_arr[i][1]-k, c_arr, n)
        if si ≤ i-2
            push!(ranges, (si,i))
        end
    end
    return ranges
end

function get_valid_ranges(N, cartesian_array, range_ind_arr, k)
    N_valid_ranges = Vector{Vector{Tuple{Int,Int}}}(undef, N);
    @inbounds for (ind,(c_start, c_end)) in enumerate(eachrow(range_ind_arr))
        c_arr = @view cartesian_array[c_start:c_end];
        N_valid_ranges[ind] = get_valid_Triplet_ranges(c_arr, k);
    end
    return N_valid_ranges
end

function record_tuple!(H, N_valid_ranges, n, c, range_ind_arr)
    @inbounds for vtr in N_valid_ranges[n]
        for i = vtr[1]:vtr[2]-2
            for j = i+1:vtr[2]-1
                c1,c2,c3 = (i,j,vtr[2]) .+ (n==1 ? 0 : range_ind_arr[n-1,2])
                m1,m2,m3,d12,d13 = c[c1][2], c[c2][2], c[c3][2], c[c2][1]-c[c1][1], c[c3][1]-c[c1][1]
                H[m1,m2,m3,d12+1,d13+1] += 1
            end
        end
    end
end

function make_enriched_hypercube(N_valid_ranges, cartesian_array, M, N, range_ind_arr, k)
    # TODO: probably a more efficient data structure exist; think about it later on
    H = zeros(UInt16, (M,M,M,k+1,k+1));
    for n = 1:N
        record_tuple!(H, N_valid_ranges, n, cartesian_array, range_ind_arr) 
    end
    return H
end

function get_number_of_outliers_via_esd(H; 
                                        count_thresh=10,
                                        r=1000, 
                                        alpha=0.01)
    #= Find all triplets whose count is more than count-thresh, perform generalized 
       ESD test to determine the number of outliers, use this number to determine 
       the count-threshold. Returns the cartesian indices in enriched hypercube H
       (each of those index has counts ≥ count-threshold)
    =#
    all_counts_ind = findall(H .> count_thresh)
    counts = Array{eltype(H),1}(undef, length(all_counts_ind));
    @inbounds for (ind,i) in enumerate(all_counts_ind)
        counts[ind] = H[i];
    end

    c̄, std_c            = count_mean_std(counts)
    Ris                 = calculate_Rᵢ(counts, c̄, std_c, r)
    lambdais            = calculate_λᵢ(counts, alpha, r)
    num_outliers        = get_number_of_outliers(Ris, lambdais)
    num_outliers == 0 && return nothing
    sort_inds           = sortperm(counts, rev=true)
    esd_count_thresh    = @view counts[sort_inds[1:num_outliers]][end]
    found = findall(H .≥ esd_count_thresh)
    return found
end

function add_key!(d, c)
    @inbounds haskey(d[c[3]], c[2]) ? push!(d[c[3]][c[2]], c[1]) : d[c[3]][c[2]]=[c[1]];
end

function get_cartesian_dict(cartesian_array, N)
    d=[Dict{Int, Vector{Int}}() for _ = 1:N];
    for c in cartesian_array
        add_key!(d, c)
    end
    return d
end

function pos_incr!(pos, ind, flen, m1, m2, m3, d12, d13, n, d)
    @inbounds if haskey(d[n], m1) && haskey(d[n], m2) && haskey(d[n], m3)
       for m1_ind in d[n][m1]
            for m2_ind in d[n][m2]                
                for m3_ind in d[n][m3]
                    (m2_ind - m1_ind == d12 && m3_ind - m1_ind == d13) && (
                        pos[ind] = (m1_ind:m3_ind+flen-1, n, false);
                        #= the third component indicate use_comp, for which 
                           we don't use but keep it here for the triplets. 
                        =#
                    ind += 1)
                end
            end
        end
    end
    return ind
end

function get_this_pos(H, m1, m2, m3, d12, d13, flen, d, N)
    pos = vec_tup_t_(undef, H[m1,m2,m3,d12,d13]);    
    ind = 1; d12_ = d12-1; d13_ = d13-1;
    for n=1:N
        ind = pos_incr!(pos, ind, flen, m1, m2, m3, d12_, d13_, n, d)
    end
    return pos
end

function obtain_positions(H, found, flen, d, N)
    new_pos = Vector{vec_tup_t_}(undef, length(found));
    @inbounds for (ind,c) in enumerate(found)
        new_pos[ind] = get_this_pos(H, c[1], c[2], c[3], c[4], c[5], flen, d, N);
    end
    return new_pos
end

function obtain_pos(Z_nz, thresh, flen; k=15, count_thresh=10, r=1000, alpha=0.01)
    _, M, N = size(Z_nz)
    # extract the highly activated sparse code components
    cartesian_array, range_ind_arr = get_indicator_arrays(Z_nz, thresh);
    # get the valid ranges for triples using binary search
    N_valid_ranges = get_valid_ranges(N, cartesian_array, range_ind_arr, k);
    # make the enriched hypercube
    H = make_enriched_hypercube(N_valid_ranges, cartesian_array, M, N, range_ind_arr, k);
    # get the cartesian indices (m1,m2,m3,d12,d13) for the highly enriched triplets 
    found = get_number_of_outliers_via_esd(H; count_thresh = count_thresh, r = r, alpha = alpha);
    isnothing(found) && return nothing
    # make the dictionary on sparse code components for fast look up
    d = get_cartesian_dict(cartesian_array, N);
    # obtain the positions for each enrichlet triplets
    new_pos = obtain_positions(H, found, flen, d, N);
    return new_pos
end


