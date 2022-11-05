"""
data: (to be added)
pos: positions and the sequence; a vector of 
    [(start₁:end₁, seq_num₁), (start₂, end₂, seq_num₂),...]
ps: pseudocount

return a position frequency matrix, i.e.,
       1   2   3   4
    A 0.1 0.1 .......
    C 0.2 0.1 .......
    G 0.3 0.1 .......
    T 0.4 0.7 .......
"""
function str2freq(data, pos, k; ps=0.01)
    count_mat = countmat(get_strs(data, pos, k)) .+ ps;
    count_mat ./ sum(count_mat, dims=1);
end

function str2count2ic(data, pos, k; ps=0.01)
   freq_mat = str2freq(data, pos, k; ps=ps)
   ic_cols = fill(0.0, size(freq_mat, 2));
   for i = 1:size(freq_mat, 2)
       ic_cols[i] = 2 + sum((@view freq_mat[:,i]) .* log2.(@view freq_mat[:,i]))
   end
   return ic_cols
end

function countvec2ic(count_vec, pos_ki; ps=0.01)  
    count_vec = count_vec .+ ps;
    count_vec_sum = sum(count_vec);
    freq_vec = count_vec ./ count_vec_sum
    count_vec_sum/length(pos_ki), 2 + sum(freq_vec .* log2.(freq_vec))
end

get_expand_ind_left(p, dec)  = p[1][1]-dec;
get_expand_ind_right(p, inc) = p[1][end]+inc;
get_shrink_ind_left(p, inc)  = p[1][1]+inc;
get_shrink_ind_right(p, dec) = p[1][end]-dec;

# count_vec_at_pos(p, data, char_ind) = char_ind < 1 || char_ind > data.L ? 
#     atcg2dummy['z'] : atcg2dummy[data.raw_data[p[2]].str[char_ind]]

count_vec_at_pos(p, data, char_ind) = char_ind < 1 || char_ind > data.L ? 
    atcg2dummy['z'] : atcg2dummy[data.raw_data[p[2]][char_ind]]

function left_expansion(data, pos_ki, dec; ic_expand_thresh=0.5, pc=1.0, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos(p, data, get_expand_ind_left(p, dec))
    end
    percentage_used, ic = countvec2ic(count_vec .+ pc, pos_ki)
    return percentage_used > tol && ic > ic_expand_thresh
end

function right_expansion(data, pos_ki, inc; ic_expand_thresh=0.5, pc=1.0, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos(p, data, get_expand_ind_right(p, inc))
    end
    percentage_used, ic = countvec2ic(count_vec .+ pc, pos_ki)
    return percentage_used > tol && ic > ic_expand_thresh
end

function left_shrinkage(data, pos_ki, inc; ic_shrink_thresh=0.2, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos(p, data, get_shrink_ind_left(p, inc))    
    end
    percentage_used, ic = countvec2ic(count_vec, pos_ki)
    # println("ic: $ic")
    return percentage_used > tol && ic < ic_shrink_thresh
end

function right_shrinkage(data, pos_ki, dec; ic_shrink_thresh=0.2, tol=0.975)
    count_vec = @SVector zeros(Float64, 4)
    for p in pos_ki
        count_vec += 
            count_vec_at_pos(p, data, get_shrink_ind_right(p, dec))    
    end
    percentage_used, ic = countvec2ic(count_vec, pos_ki)
    # isnan(ic) && println(count_vec)
    return percentage_used > tol && ic < ic_shrink_thresh
end

function expansion_left_right(data, pos_ki, ic_expand_thresh)
    expand_left = true; expand_right = true; 
    left_dec = 1; right_inc = 1;
    while expand_left 
        expand_left = 
            left_expansion(data, pos_ki, left_dec; ic_expand_thresh=ic_expand_thresh)
        expand_left ? left_dec+=1 : left_dec-=1;
    end
    while expand_right
        expand_right = 
            right_expansion(data, pos_ki, right_inc; ic_expand_thresh=ic_expand_thresh)
        expand_right ? right_inc+=1 : right_inc-=1;
    end
    return left_dec, right_inc
end

# TODO: can acutally form count matrices and trim accordingly

function trimming_left_right(expansion::Tuple{Int, Int}, data, pos_ki, ic_shrink_thresh)
    # returns how much increment (decrement) from the left (from the right)
    shrink_left_  = expansion[1] == 0 ? true : false; 
    shrink_right_ = expansion[2] == 0 ? true : false;
    left_inc = shrink_left_ ? 1 : 0; 
    right_dec = shrink_right_ ? 1 : 0;
    while shrink_left_
        shrink_left_ = 
            left_shrinkage(data, pos_ki, left_inc-1; ic_shrink_thresh=ic_shrink_thresh)
        shrink_left_ ? left_inc+=1 : left_inc -= 1;
    end
    while shrink_right_
        shrink_right_ = 
            right_shrinkage(data, pos_ki, right_dec-1; ic_shrink_thresh=ic_shrink_thresh)
        shrink_right_ ? right_dec+=1 : right_dec-=1;
    end
    return left_inc, right_dec
end

function msa_expansion(pos, data, ic_expand_thresh, ic_shrink_thresh)
    expansions = Vector{Tuple{Int,Int}}();
    shrinkage = Vector{Tuple{Int,Int}}();
    for pos_ki in pos
        push!(expansions, expansion_left_right(data, pos_ki, ic_expand_thresh))
    end
    for (ind, expansion) in enumerate(expansions)
        push!(shrinkage, trimming_left_right(expansion, data, pos[ind], ic_shrink_thresh))
    end
    return expansions, shrinkage
end

function get_new_range(p_i, el, er, sl, sr)
    (p_i[1][1]-el+sl):(p_i[1][end]+er-sr), p_i[2], false
end

function get_new_range!(v, p_i, L, el, er, sl, sr)
    new_range = get_new_range(p_i, el, er, sl, sr)
    (new_range[1][1] ≥ 1 && new_range[1][end] ≤ L) && push!(v, new_range)
end

function pos_range_edit(expand__left, expand__right, 
                        shrink__left, shrink__right, pos_ki, L)
    @assert expand__left == 0 || shrink__left == 0 "one left shink/expand has to be 0"
    @assert expand__right == 0 || shrink__right == 0 "one right shink/expand has to be 0"    

    v = vec_tup_t_();
    for p in pos_ki
        get_new_range!(v, p, L, expand__left, expand__right, shrink__left, shrink__right);
    end
    return v
end

function get_expanded_pos(pos, data; 
                          ic_expand_thresh=0.5, ic_shrink_thresh=0.4, min_len=4)
    expansions, shrinkage = 
        msa_expansion(pos, data, ic_expand_thresh, ic_shrink_thresh)
    new_pos = Vector{eltype(pos)}();
    for (ind, ((el,er),(sl,sr))) in enumerate(zip(expansions, shrinkage))
        # if (get_pos_len(pos[ind])-sr+1)-sl > 0            
        if get_pos_len(pos[ind])-sr-sl > min_len
            push!(new_pos, pos_range_edit(el, er, sl, sr, pos[ind], data.L))
        # else
            # @info "len: $(get_pos_len(pos[ind])-sr-sl) too small"
            # @info "len: $((get_pos_len(pos[ind])-sr+1)-sl) too small"
        end
    end
    return new_pos
end