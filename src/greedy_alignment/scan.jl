function motifs_prep(ms::motifs{T,S}) where {T<:Integer,S<:Real}
    positions = [Dict{T,Vector{T}}()    for _ = 1:ms.num_motifs];
    scores    = [Dict{T,Vector{S}}()    for _ = 1:ms.num_motifs];
    use_comp  = [Dict{T,Vector{Bool}}() for _ = 1:ms.num_motifs];
    return positions, scores, use_comp
end

function scan_n!(n, data_n, pwm, pwm_comp, pwm_len, L, p_k, s_k, u_k, thresh_k)
    # scan the sequence in the orientation as in the dataset 
    @inbounds for i = 1:L-pwm_len+1 # devectorization
        score = eltype(pwm)(0); score_c = eltype(pwm_comp)(0);
        for j = 1:4, k = 1:pwm_len
            score += pwm[j,k] * data_n[j,i+k-1];
            score_c += pwm_comp[j,k] * data_n[j,i+k-1];
        end
        if score > thresh_k 
            if haskey(p_k, n)
                push!(p_k[n], i); push!(s_k[n], score); push!(u_k[n], false)
            else
                p_k[n] = [i]; s_k[n] = [score]; u_k[n] = [false];
            end
        end
        if score_c > thresh_k
            if haskey(p_k, n)
                push!(p_k[n], i); push!(s_k[n], score_c); push!(u_k[n], true)
            else
                p_k[n] = [i]; s_k[n] = [score_c]; u_k[n] = [true];
            end
        end 
    end
end

function overlapping_scan!(ms::motifs, data; test=false)
    data_matrix = test ? data.data_matrix_test : data.data_matrix;
    N = test ? data.N_test : data.N;

    p, s, u = motifs_prep(ms);

    @inbounds for k = 1:ms.num_motifs   
        pwm_comp_k = mat_complement(ms.pwms[k]); 
        for n = 1:N
            data_n = reshape(Base.view(data_matrix,:,1,n), (4, data.L));
            scan_n!(n, data_n, ms.pwms[k], pwm_comp_k, 
                ms.lens[k], data.L, p[k], s[k], u[k], ms.thresh[k]);
        end
    end
  
    ms.positions = p;
    ms.scores = s;
    ms.use_comp = u;
end

function overlapping_scan_both!(ms::motifs, data; test=false)
    data_matrix = cat(data.data_matrix_test, data.data_matrix, dims=3)
    N = data.N+data.N_test;

    p, s, u = motifs_prep(ms);

    @inbounds for k = 1:ms.num_motifs   
        pwm_comp_k = mat_complement(ms.pwms[k]); 
        for n = 1:N
            data_n = reshape(Base.view(data_matrix,:,1,n), (4, data.L));
            scan_n!(n, data_n, ms.pwms[k], pwm_comp_k, 
                ms.lens[k], data.L, p[k], s[k], u[k], ms.thresh[k]);
        end
    end
  
    ms.positions = p;
    ms.scores = s;
    ms.use_comp = u;
end

function overlapping_scan_bg!(ms::motifs, data; test=false)
    data_matrix_bg = test ? data.data_matrix_bg_test : data.data_matrix_bg;

    p, s, u = motifs_prep(ms);
    @inbounds for k = 1:ms.num_motifs   
        pwm_comp_k = mat_complement(ms.pwms[k]);      
        for n = 1:size(data_matrix_bg,2)
            data_n = reshape(Base.view(data_matrix_bg,:,n), (4, data.L)); 
            scan_n!(n, data_n, ms.pwms[k], pwm_comp_k, 
                ms.lens[k], data.L, p[k], s[k], u[k], ms.thresh[k]);
        end
    end
    ms.positions_bg = p;
    ms.scores_bg = s;
    ms.use_comp_bg = u;
end

function max_score_ind(scores_n, scores_n_mask)
    max_score = -Inf;
    maxscore_ind = nothing;
    maxscore_ind_m = nothing;
    @inbounds for (ind_1,s) in enumerate(scores_n)
        for (ind_2, s_) in enumerate(s)
            if s_ > max_score && scores_n_mask[ind_1][ind_2]
                max_score = s_
                maxscore_ind = ind_2; 
                maxscore_ind_m = ind_1;
            end
        end
    end
    return maxscore_ind, maxscore_ind_m
end

max_pos_ind(positions_n, max_score_ind, max_score_m) = 
    positions_n[max_score_m][max_score_ind]

function non_overlap_scan!(ms::motifs{T,S}, data; test=false) where {T <: Int, S <: Real}
    N = test ? data.N_test : data.N;    
    @inbounds for n = 1:N
        # println(n)
        spans_pos = T[]; 
        spans_len = T[];

        indices_to_keep_n = [
            haskey(ms.positions[i], n) ? fill(false, length(ms.positions[i][n])) : Bool[]
             for i = 1:ms.num_motifs];
        scores_n  = [haskey(ms.scores[i], n) ? ms.scores[i][n] : S[] for i = 1:ms.num_motifs];
        scores_n_mask  = [fill(true, length(s)) for s in scores_n]; # so entries in ms.scores aren't modified
        positions_n = [haskey(ms.positions[i], n) ? ms.positions[i][n] : T[] for i = 1:ms.num_motifs];   
        maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)

        if !isnothing(maxscore_ind)            
            while !isnothing(maxscore_ind)
                maxpos_ind = max_pos_ind(positions_n, maxscore_ind, maxscore_ind_m)
                intersect_ = false;
                for (p,l) in zip(spans_pos, spans_len)
                    # check whether this "segment" intersect with any previous segments
                    p_end     = p+l-1;
                    p_max     = positions_n[maxscore_ind_m][maxscore_ind];
                    p_max_end = p_max+ms.lens[maxscore_ind_m]-1;
                    if p ≤ p_max ≤ p_end || 
                        p ≤ p_max_end ≤ p_end || 
                        p_max ≤ p ≤ p_max_end || 
                        p_max ≤ p_end ≤ p_max_end
                        intersect_ = true;
                        break
                    end
                end
                if !intersect_
                    indices_to_keep_n[maxscore_ind_m][maxscore_ind] = true;
                    push!(spans_pos, maxpos_ind);
                    push!(spans_len, ms.lens[maxscore_ind_m]);
                end
                scores_n_mask[maxscore_ind_m][maxscore_ind]=false;
                maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)
            end
            for j = 1:ms.num_motifs         
                if haskey(ms.positions[j], n)
                    ms.positions[j][n] = ms.positions[j][n][indices_to_keep_n[j]];
                    ms.scores[j][n] = ms.scores[j][n][indices_to_keep_n[j]];
                    ms.use_comp[j][n] = ms.use_comp[j][n][indices_to_keep_n[j]];
                    # so that the set of sequences covered by pwms are disjoint
                end
            end
        end
    end
end

function non_overlap_scan_both!(ms::motifs{T,S}, data; test=false) where {T <: Int, S <: Real}
    N = data.N+data.N_test;
    @inbounds for n = 1:N
        # println(n)
        spans_pos = T[]; 
        spans_len = T[];

        indices_to_keep_n = [
            haskey(ms.positions[i], n) ? fill(false, length(ms.positions[i][n])) : Bool[]
             for i = 1:ms.num_motifs];
        scores_n  = [haskey(ms.scores[i], n) ? ms.scores[i][n] : S[] for i = 1:ms.num_motifs];
        scores_n_mask  = [fill(true, length(s)) for s in scores_n]; # so entries in ms.scores aren't modified
        positions_n = [haskey(ms.positions[i], n) ? ms.positions[i][n] : T[] for i = 1:ms.num_motifs];   
        maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)

        if !isnothing(maxscore_ind)            
            while !isnothing(maxscore_ind)
                maxpos_ind = max_pos_ind(positions_n, maxscore_ind, maxscore_ind_m)
                intersect_ = false;
                for (p,l) in zip(spans_pos, spans_len)
                    # check whether this "segment" intersect with any previous segments
                    p_end     = p+l-1;
                    p_max     = positions_n[maxscore_ind_m][maxscore_ind];
                    p_max_end = p_max+ms.lens[maxscore_ind_m]-1;
                    if p ≤ p_max ≤ p_end || 
                        p ≤ p_max_end ≤ p_end || 
                        p_max ≤ p ≤ p_max_end || 
                        p_max ≤ p_end ≤ p_max_end
                        intersect_ = true;
                        break
                    end
                end
                if !intersect_
                    indices_to_keep_n[maxscore_ind_m][maxscore_ind] = true;
                    push!(spans_pos, maxpos_ind);
                    push!(spans_len, ms.lens[maxscore_ind_m]);
                end
                scores_n_mask[maxscore_ind_m][maxscore_ind]=false;
                maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)
            end
            for j = 1:ms.num_motifs         
                if haskey(ms.positions[j], n)
                    ms.positions[j][n] = ms.positions[j][n][indices_to_keep_n[j]];
                    ms.scores[j][n] = ms.scores[j][n][indices_to_keep_n[j]];
                    ms.use_comp[j][n] = ms.use_comp[j][n][indices_to_keep_n[j]];
                    # so that the set of sequences covered by pwms are disjoint
                end
            end
        end
    end
end

function non_overlapping_scan_bg!(ms::motifs{T,S}, data; test=false) where {T,S}
    N = test ? data.N_test : data.N;
    @inbounds for n = 1:N
        spans_pos = T[]; 
        spans_len = T[];

        indices_to_keep_n = [
            haskey(ms.positions_bg[i], n) ? fill(false, length(ms.positions_bg[i][n])) : Bool[]
             for i = 1:ms.num_motifs];
        scores_n  = [haskey(ms.scores_bg[i], n) ? ms.scores_bg[i][n] : S[] for i = 1:ms.num_motifs];
        scores_n_mask  = [fill(true, length(s)) for s in scores_n]; # so entries in scores aren'
        positions_n = [haskey(ms.positions_bg[i], n) ? ms.positions_bg[i][n] : T[] for i = 1:ms.num_motifs];   
        maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)

        if !isnothing(maxscore_ind)     
            while !isnothing(maxscore_ind)
                maxpos_ind = max_pos_ind(positions_n, maxscore_ind, maxscore_ind_m)
                intersect_ = false;
                for (p,l) in zip(spans_pos,spans_len)
                    # check whether this "segment" intersect with any previous segments
                    p_end     = p+l-1;
                    p_max     = positions_n[maxscore_ind_m][maxscore_ind];
                    p_max_end = p_max+ms.lens[maxscore_ind_m]-1;
                    if p ≤ p_max ≤ p_end || 
                        p ≤ p_max_end ≤ p_end || 
                        p_max ≤ p ≤ p_max_end || 
                        p_max ≤ p_end ≤ p_max_end
                        intersect_ = true;
                    end
                end
                if !intersect_
                    indices_to_keep_n[maxscore_ind_m][maxscore_ind] = true;
                    push!(spans_pos, maxpos_ind);
                    push!(spans_len, ms.lens[maxscore_ind_m]);                    
                end
                scores_n_mask[maxscore_ind_m][maxscore_ind]=false;
                maxscore_ind, maxscore_ind_m = max_score_ind(scores_n, scores_n_mask)
            end
            for j = 1:ms.num_motifs         
                if haskey(ms.positions_bg[j], n)
                    ms.positions_bg[j][n] = ms.positions_bg[j][n][indices_to_keep_n[j]];
                    ms.scores_bg[j][n] = ms.scores_bg[j][n][indices_to_keep_n[j]];
                    # so that the set of sequences covered by pwms are disjoint
                end
            end
        end
    end
end