function get_activate_counts_p(p)
    sum_p = 0
    for key in keys(p)
        for r in p[key]
            sum_p += r[end]-r[1]+1
        end
    end
    return sum_p
end

function get_activate_counts(ms::motifs; bg = false)
    positions = bg ? ms.positions_bg : ms.positions;
    activate_count = Vector{Float64}(undef, length(ms.lens))
    for i = 1:ms.num_motifs
        sum_p = 0
        for key in keys(positions[i])
            sum_p += length(positions[i][key])
        end
        activate_count[i] = sum_p
        # println("Motif $i: $sum_p")
    end
    return activate_count
end

function get_activate_counts_use_comp(ms::motifs; comp=true)
    activate_count = Vector{Float64}(undef, length(ms.lens))
    for i = 1:ms.num_motifs
        sum_p = 0
        for key in keys(ms.use_comp[i])
            sum_p += sum(ms.use_comp[i][key] .== comp)
        end
        activate_count[i] = sum_p
        # println("Motif $i: $sum_p")
    end
    return activate_count
end

function get_activate_counts_bg(ms::motifs)
    activate_count_bg = Vector{Float64}(undef, length(ms.lens))
    for i = 1:ms.num_motifs
        sum_p = 0
        for key in keys(positions_bg[i])
            sum_p += length(positions_bg[i][key])
        end
        activate_count_bg[i] = sum_p
        # println("Motif $i: $sum_p")
    end
    return activate_count_bg
end

# get the activate counts in the control data set
get_both_activate_counts(ms::motifs) =
    get_activate_counts(ms),  get_activate_counts(ms; bg=true)

function get_fisher_p_values(ms::motifs, data; test=false)
    (!test && isnothing(ms.positions)) && scan_trainset_cpu!(ms, data, true);
    test && scan_testset!(ms, data);
    activate_counts, activate_counts_bg = get_both_activate_counts(ms);
    activate_sum = (test ? data.N_test : data.N) * data.L # total number of components that can be activated
    pvalues = fill(0.0, ms.num_motifs);  

    for i = 1:ms.num_motifs
        a = activate_counts[i]; b = activate_counts_bg[i];
        c = activate_sum - a; d = activate_sum - b;
        if a == 0 && b == 0
            pvalues[i] = 1.0 # if there's no activation in both just throw it away
        else
            q = FisherExactTest(promote_i(a, c, b, d)...);
            pvalues[i] = HypothesisTests.pvalue(q, tail=:right);
        end
    end
    return pvalues
end