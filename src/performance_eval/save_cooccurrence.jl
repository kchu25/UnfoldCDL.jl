function fill_LMN_cube(data, ms, N)
    # fill the cube with activation indicators
    LMN_cube_pos  = falses(data.L, ms.num_motifs, N);
    LMN_cube_comp = falses(data.L, ms.num_motifs, N);
    # fill
    @inbounds for m = 1:ms.num_motifs
        if !isnothing(ms.positions[m])
            for k in keys(ms.positions[m])
                for (ind,pos) in enumerate(ms.positions[m][k])
                    LMN_cube_pos[pos,m,k] = 1; 
                    LMN_cube_comp[pos,m,k] = ms.use_comp[m][k][ind] ? 1 : 0;
                end
            end
        end
    end
    return LMN_cube_pos, LMN_cube_comp
end

function add_pos_diff!(m_front, m_back, 
                       pos_front, pos_back, 
                       comp_front, comp_back, 
                       dict, len_front, n, a, alpha)
    d = pos_back - (pos_front+len_front-1);
    k = (m_front, m_back, comp_front, comp_back);
    if d ≤ alpha
        @inbounds if haskey(dict, k)
            if haskey(dict[k], (n,a))
                push!(dict[k][n,a], d);
            else
                dict[k][n,a] = [d];
            end
        else
            dict[k] = Dict((n,a)=>[d]);
        end
    end
end

function fill_dict_inner_inner!(ms, dict, active_pos_m1n, active_pos_m2n, 
                                use_comp_m1, use_comp_m2, m1, m2, n, alpha)
    @inbounds for (ind1,a1) in enumerate(active_pos_m1n), (ind2,a2) in enumerate(active_pos_m2n)
        if a1 > a2
            add_pos_diff!(m2, m1, a2, a1, 
                            use_comp_m2[ind2], use_comp_m1[ind1], dict, ms.lens[m2], n, a1, alpha)
        elseif a2 > a1
            add_pos_diff!(m1, m2, a1, a2, 
                            use_comp_m1[ind1], use_comp_m2[ind2], dict, ms.lens[m1], n, a2, alpha)
        end
    end
end

function return_active_pos(m1, m2, n, LMN_cube_pos, LMN_cube_comp)
    active_pos_m1n = findall((@view LMN_cube_pos[:,m1,n]) .== 1);
    use_comp_m1    = @view LMN_cube_comp[:,m1,n][active_pos_m1n];
    @inbounds if m1 == m2 
        return active_pos_m1n, active_pos_m1n, use_comp_m1, use_comp_m1
    else
        active_pos_m2n = findall((@view LMN_cube_pos[:,m2,n]) .== 1);
        use_comp_m2    = @view LMN_cube_comp[:,m2,n][active_pos_m2n];
        return active_pos_m1n, active_pos_m2n, use_comp_m1, use_comp_m2
    end
end

function fill_dict_inner!(data, ms, m1, m2, LMN_cube_pos, LMN_cube_comp, dict, alpha)
    for n = 1:(data.N+data.N_test)
        # println("m1: $m1, m2: $m2, n: $n")
        active_pos_m1n, active_pos_m2n, use_comp_m1, use_comp_m2 = 
            return_active_pos(m1, m2, n, LMN_cube_pos, LMN_cube_comp)
        fill_dict_inner_inner!(ms, dict, active_pos_m1n, active_pos_m2n, 
            use_comp_m1, use_comp_m2, m1, m2, n, alpha)
    end
end

function fill_dict(ms, data, LMN_cube_pos, LMN_cube_comp, alpha)
    dict = Dict{Tuple{Int,Int,Bool,Bool}, Dict{Tuple{Int,Int}, Vector{Int}}}();
    # dict = Dict{Tuple{Int,Int,Bool,Bool}, Dict{Tuple{Int,Int}, Int}}();
    @inbounds for m1 = 1:ms.num_motifs, m2 = m1:ms.num_motifs
        fill_dict_inner!(data, ms, m1, m2, LMN_cube_pos, LMN_cube_comp, dict, alpha)
    end
    return dict
end

function prune_dict(dict)
    new_dict = Dict{Tuple{Int,Int,Bool,Bool}, Vector{Int}}();
    for k in keys(dict)
        new_dict[k] = Int[]
        for q in keys(dict[k])
            # push!(new_dict[k], dict[k][q])
            append!(new_dict[k], dict[k][q])
        end
    end
    new_dict
end

function get_merged_d_arr(c, cc, dict)
    haskey_c = haskey(dict, c);
    haskey_cc = haskey(dict, cc);
    (haskey_c && haskey_cc) && return append!(dict[c], dict[cc]);
    (haskey_c && !haskey_cc) && return dict[c];
    (!haskey_c && haskey_cc) && return dict[cc];
    (!haskey_c && !haskey_cc) && return Int[];
end

function fill_dict_merged(ms, dict)
    # dict_merged = Dict{Tuple{Int,Int,Bool,Bool}, Vector{Tuple{Int,Int}}}();
    dict_merged = Dict{Tuple{Int,Int,Bool,Bool}, Vector{Int}}();
    @inbounds for m1 = 1:ms.num_motifs, m2 = m1:ms.num_motifs        
        c1 = (m1,m2,true,true);     c1c = (m2,m1,false,false);
        c2 = (m1,m2,true,false);    c2c = (m2,m1,true,false);
        c3 = (m1,m2,false,true);    c3c = (m2,m1,false,true);
        c4 = (m1,m2,false,false);   c4c = (m2,m1,true,true);    
        if m1 == m2
            dict_merged[c1] = get_merged_d_arr(c1,c1c,dict);
            dict_merged[c2] = get_merged_d_arr(c2,c2c,dict);
            dict_merged[c3] = get_merged_d_arr(c3,c3c,dict);
        else
            dict_merged[c1] = get_merged_d_arr(c1,c1c,dict);
            dict_merged[c2] = get_merged_d_arr(c2,c2c,dict);
            dict_merged[c3] = get_merged_d_arr(c3,c3c,dict);
            dict_merged[c4] = get_merged_d_arr(c4,c4c,dict);
        end
    end
    return dict_merged
end

function get_pmi_counts(dict_merged, ms)
    normalize_factor = 0;
    marginal_count_left  = zeros(ms.num_motifs);
    marginal_count_right = zeros(ms.num_motifs);
    @inbounds for k in keys(dict_merged)
        k1, k2, _, _ = k;
        marginal_count_left[k1] += length(dict_merged[k]);
        marginal_count_right[k2] += length(dict_merged[k]);
        normalize_factor += length(dict_merged[k])
    end
    return normalize_factor, marginal_count_left, marginal_count_right
end

function get_pmi(dict_merged, 
                 normalize_factor, 
                 marginal_count_left, 
                 marginal_count_right
                 )
    pmis = Dict(k=>0.0 for k in keys(dict_merged));
    @inbounds for k in keys(dict_merged)
        k1, k2, _, _ = k;
        pmis[k] = log2(
            (length(dict_merged[k]) * normalize_factor) / 
                (marginal_count_left[k1]*marginal_count_right[k2]) )
    end
    return pmis
end

function get_config_used_vec(pmis, use_vec)
    pmis_keys2condsider = [k for k in keys(pmis) if  pmis[k] > 0];
    config_used_vec = [i for i in pmis_keys2condsider if use_vec[i[1]] && use_vec[i[2]]];
    return config_used_vec
end

function filter_dict_merge_w_alpha!(dict_merged, alpha)
    for k in keys(dict_merged)
        dict_merged[k] = dict_merged[k][findall(dict_merged[k] .≤ alpha)]
    end
end

function remap_filter_num(sort_perm, config_used_vec)
    map_ms_motif_num = Dict(i=>ind for (ind, i) in enumerate(sort_perm))
    return [(map_ms_motif_num[m1], map_ms_motif_num[m2], u1, u2) 
                for (m1,m2,u1,u2) in config_used_vec]
end

function plot_cooccurrence_png(dict_merged, remapped_key,
                               config_used_vec, npmis, co_occ_count,
                               pwms_img, pwms_img_c,
                               alpha,
                               save_where; 
                               text_size=24,
                               font="TeX Gyre Heros Bold",
                               colormap = ColorSchemes.glasbey_hv_n256.colors,
                               with_words=true
                               )
                               
    K = length(config_used_vec);
    yticks = [1:10:maximum(dict_merged[config_used_vec[k]]) for k = 1:K];
    f = Figure(resolution=(1700, K*250));
    my_grid = f[1:K,1:3] = GridLayout();
    plot_axs = [Axis(my_grid[k,2], xticklabelsize=22, xlabelsize=32, width=850, 
                        xlabelfont=font, yticks=yticks[k]) for k = 1:K];
    left_motif_axs  = [Axis(my_grid[k,1], 
               spinewidth=0, aspect=DataAspect(), xreversed=true, 
               xgridvisible=false, ygridvisible=false,
               leftspinevisible=false, bottomspinevisible=false, 
               rightspinevisible=false, topspinevisible=false, 
               xticksvisible=false, yticksvisible=false, 
               xticklabelsvisible=false, yticklabelsvisible=false) 
               for k=1:K];
    right_motif_axs = [Axis(my_grid[k,3],
               spinewidth=0, aspect=DataAspect(), halign=:right, 
               leftspinevisible=false, bottomspinevisible=false,  
               xgridvisible=false, ygridvisible=false,
               rightspinevisible=false, topspinevisible=false, 
               xticksvisible=false, yticksvisible=false, 
               xticklabelsvisible=false, yticklabelsvisible=false)
               for k=1:K];

    @inbounds for k = 1:K
        k < K && linkxaxes!(right_motif_axs[k], right_motif_axs[k+1])
        k < K && linkxaxes!(left_motif_axs[k], left_motif_axs[k+1])
        m1, m2, uc1, uc2 = remapped_key[k]
        if uc1 == 1 && uc2 == 1
            image!(left_motif_axs[k], rotr90(reverse(pwms_img[m2], dims=2)))
            image!(right_motif_axs[k], rotr90(pwms_img[m1]))  

            left_motif_axs[k].xlabel = "D$m2";
            right_motif_axs[k].xlabel = "D$m1";
        else
            image!(left_motif_axs[k],  uc1 ? rotr90(reverse(pwms_img_c[m1], dims=2)) : rotr90(reverse(pwms_img[m1],dims=2)))
            image!(right_motif_axs[k], uc2 ? rotr90(pwms_img_c[m2]) : rotr90(pwms_img[m2]));

            left_motif_axs[k].xlabel  = uc1 ? "D$m1-rev.comp" : "D$m1";
            right_motif_axs[k].xlabel = uc2 ? "D$m2-rev.comp" : "D$m2";
        end
        left_motif_axs[k].xlabelsize=text_size-1;
        right_motif_axs[k].xlabelsize=text_size-1;
    end
   
    @inbounds for k = 1:K-1 linkxaxes!(plot_axs[k], plot_axs[k+1]) end
    @inbounds for k = 1:K-1 hidedecorations!(plot_axs[k], grid=false) end

    @inbounds for k = 1:K
        category_labels = fill("", length(dict_merged[config_used_vec[k]]))
        # hist!(plot_axs[k], bins=alpha, dict_merged[config_used_vec[k]], 
        #       strokewidth = 2, strokecolor = :black)
        rainclouds!(plot_axs[k], category_labels, dict_merged[config_used_vec[k]];
            xlabel = "number of nucleotides in between",
            orientation = :horizontal, jitter_width=0.01, 
            plot_boxplots = true, cloud_width=0.5, clouds=violin,
            markersize=0, side_nudge=0, 
            color = colormap[indexin(category_labels, unique(category_labels))])
        if with_words
            text!(plot_axs[k], "nPMI: $(npmis[k]) ", 
                position = (alpha-3, 1.17), textsize=text_size)
            text!(plot_axs[k], "# co-occ: $(co_occ_count[k]) ", 
                position = (alpha-3, 1.13), textsize=text_size)
        end
    end

    for k = 1:K
        linkyaxes!(left_motif_axs[k], right_motif_axs[k]);
    end
    if with_words
        save(save_where*"/coocc_$(alpha)_ww.png", f, px_per_unit = 1)
    else
        save(save_where*"/coocc_$alpha.png", f, px_per_unit = 1)
    end
end

function plot_cooccurrence(ms, 
                           data, 
                           sort_perm, 
                           use_vec, 
                           pwms_img, 
                           pwms_img_c, 
                           save_loc::String; 
            alphas=[-1, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
            )
    LMN_cube_pos, LMN_cube_comp = 
        fill_LMN_cube(data, ms, data.N+data.N_test);
    alpha_ind = length(alphas); alpha = alphas[alpha_ind];
    dict = fill_dict(ms, data, LMN_cube_pos, LMN_cube_comp, alpha);
    dict_pruned = prune_dict(dict);
    dict_merged = fill_dict_merged(ms, dict_pruned);
    
    valid_alphas = Int[];

    while alpha > -1
        normalize_factor, marginal_count_left, marginal_count_right = 
            get_pmi_counts(dict_merged, ms);
        pmis = get_pmi(dict_merged, normalize_factor, marginal_count_left, marginal_count_right);
        config_used_vec = get_config_used_vec(pmis, use_vec);
        remapped_key = remap_filter_num(sort_perm, config_used_vec);
        npmis = [round(pmis[c]/(-log2(length(dict_merged[c])/normalize_factor)), digits=2) 
            for c in config_used_vec];
        co_occ_count = [length(dict_merged[c]) for c in config_used_vec];
        if length(co_occ_count) != 0
            plot_cooccurrence_png(dict_merged, remapped_key,
                          config_used_vec, npmis, co_occ_count, 
                          pwms_img, pwms_img_c,
                          alpha,
                          save_loc
                          )
            plot_cooccurrence_png(dict_merged, remapped_key,
                          config_used_vec, npmis, co_occ_count, 
                          pwms_img, pwms_img_c,
                          alpha,
                          save_loc; with_words=false
                          )
            pushfirst!(valid_alphas, alpha);
        end
        alpha_ind -= 1; alpha = alphas[alpha_ind];
        filter_dict_merge_w_alpha!(dict_merged, alpha);
    end

    return valid_alphas
end
