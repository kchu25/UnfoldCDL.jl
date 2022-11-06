function train_and_get_perf(data; 
                            M=160, flen=8, num_epochs=25, K_c=3,
                            batch_size=48, gamma_1=0.005,
                            esd_alpha=0.05, esd_r=1000, 
                            k=8, esd_count_thresh=10, percentage_thresh=0.0055,
                            ic_expand_t=0.7, ic_shrink_t=0.4, 
                            allr_thresh=0.225,
                            alpha_fisher=1e-5,
                            indep_run=8,
                            dep_run=5,
                            num_trials=5, diff_tol=4,
                            learning_rate=0.0003
                            )
    ms = nothing;
    while num_trials > 0
        try
            @info "Training uCDL..."
            Z, _ = train_basic(data.data_matrix; 
                                M=M, 
                                num_epochs=num_epochs, 
                                f=flen, 
                                K_c=K_c, 
                                gamma_1=gamma_1, 
                                batch_size=batch_size,
                                learning_rate=learning_rate);

            @info "Obtaining enriched triplets as PWMs..."
            new_cmats = obtain_count_matrices(Z, 
                                data, flen;
                                esd_alpha=esd_alpha,
                                esd_r=esd_r, k=k,
                                esd_count_thresh=esd_count_thresh,
                                percentage_thresh=percentage_thresh,
                                ic_expand_t=ic_expand_t,
                                ic_shrink_t=ic_shrink_t,
                                allr_thresh=allr_thresh,
                                diff_tol=diff_tol);
        
            @info "Greedy alignment..."
            ms = greedy_alignment(new_cmats, data;
                                        allr_thresh=allr_thresh,
                                        alpha_fisher=alpha_fisher,
                                        ic_expand_t=ic_expand_t,
                                        ic_shrink_t=ic_shrink_t,
                                        indep_run=indep_run,
                                        dep_run=dep_run,
                                        diff_tol=diff_tol,
                                        gpu=true
                                        );

            !isnothing(ms) && !(length(ms.pfms) == 0) && break
        catch y
            if isa(y, ArgumentError)
                @info "Caught Argument Error; rerun remaining times: $num_trials"
                num_trials-=1;
            else
                @info "Caught unknown error. rerun remaining times: $num_trials"
                num_trials-=1;                
            end
            M-=10; learning_rate-=0.000025; percentage_thresh+=0.0005;
            num_trials == 0 && return nothing;
        end
    end
    return ms
end

function set_output_folder(input_fasta, output_folder)
    fasta_name = split(input_fasta, "/")[end];
    output_folder = 
        output_folder[end] == '/' ? 
            output_folder*fasta_name : 
                output_folder*"/"*fasta_name;
    return output_folder*"_result"
end

function find_motif(input_fasta::String, output_folder::String; 
                    p_value_cutoff=1e-6, co_occurrence_results=true)
    output_folder = set_output_folder(input_fasta, output_folder);    
    data = FASTA_DNA{Float32}(input_fasta); 
    ms = train_and_get_perf(data);
    @info "Saving the results..."
    save_result(ms, data, output_folder, p_value_cutoff; 
        co_occurrence_results=co_occurrence_results);
end