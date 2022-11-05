############# generalized ESD test ###################################
function get_counts(H; count_thresh=10)
    all_counts_ind = findall(H .> count_thresh);
    counts = Array{eltype(H),1}(undef, length(all_counts_ind));
    @inbounds for (ind,i) in enumerate(all_counts_ind)
        counts[ind] = H[i];
    end
    return counts
end

count_mean_std(counts_) = mean(counts_), std(counts_)
get_Ri_values(counts_, c̄, std_c) = (abs.(counts_ .- c̄)) ./ std_c

function calculate_Rᵢ(counts_, c̄, std_c, r::Int)
    counts__ = copy(counts_)
    Ris = Float64[]
    # calculate R₁,…,Rᵣ
    for _ = 1:r
        abs_deviations_normalized = get_Ri_values(counts__, c̄, std_c)
        argmax_ind = argmax(abs_deviations_normalized)
        push!(Ris, abs_deviations_normalized[argmax_ind])
        deleteat!(counts__, argmax_ind)
    end
    Ris
end

function calculate_λᵢ(counts_, alpha, r::Int)
    n = length(counts_);
    lambdas = Float64[];
    for i = 1:r
        t = TDist(n-i-1); # degree of freedom ν = n-i-1
        p = 1 - alpha / (2*(n-i+1))
        t_pv = quantile(t, p);
        lambda_i = ((n-i)*t_pv) / sqrt((n-i-1+t_pv^2)*(n-i+1))
        push!(lambdas, lambda_i)
    end
    lambdas
end

get_number_of_outliers(Ris, lambdais) = sum(Ris .> lambdais)

function get_number_of_outliers(H; 
                                count_thresh=10, 
                                r=100, # upperbound for number of outliers
                                alpha=0.05     # significance level
                                )
    
    counts = get_counts(H; count_thresh=count_thresh)

    c̄, std_c = count_mean_std(counts)
    Ris = calculate_Rᵢ(counts, c̄, std_c, r)
    lambdais = calculate_λᵢ(counts, alpha, r)
    return get_number_of_outliers(Ris, lambdais)
end

######################################################################