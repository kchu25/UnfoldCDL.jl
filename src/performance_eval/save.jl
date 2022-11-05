
const logo_folder_name = "logos"
const pics_folder_name = "other_pics"
const grey_tag_front = "<p style=\"color:grey\">"
const grey_tag_back = "</p>"

function get_significant_motif_positions_pvec(ms, data, alpha_fisher)
    pval_vec = Cdlunroll.get_fisher_p_values(ms, data; test=true)
    return pval_vec, pval_vec .< alpha_fisher;
end

make_grey(s::String) = grey_tag_front*s*grey_tag_back

function get_rounded_pval(pval::Real, low_pval)
    str = "$pval"; s = nothing;
    if !occursin("e-", str)
        s =  string(round(pval, sigdigits=3));        
    else
        q = split(str, "e-"); 
        s = join([q[1][1:4], q[2]], "e-");
    end
    return !low_pval ? make_grey(s) : s;
end

function save_pfms_as_transfac(logo_folder::String, ms, sort_perm::Vector{Int})
    for (i,ind) in enumerate(sort_perm)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(ms.pfms[ind] .* 1000)); # make it a count matrix
        for j = 1:size(ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac -n 40 --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 16 --color-scheme classic -o $(logo_folder)/d$(i).png`);

        # do it for the reverse complement as well
        io = open(logo_folder*"/d$(i)_c.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(reverse(ms.pfms[ind]) .* 1000));
        for j = 1:size(ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac -n 40 --errorbars NO -F png --fineprint " " --resolution 600 -s large --fontsize 16 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end

function save_pfms_as_transfac_noXY(logo_folder::String, ms, sort_perm::Vector{Int})
    for (i,ind) in enumerate(sort_perm)
        io = open(logo_folder*"/d$i.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(ms.pfms[ind] .* 1000)); # make it a count matrix
        for j = 1:size(ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac -n 40 --errorbars NO -F png -X NO -Y NO --fineprint " " --resolution 600 -s large --fontsize 32 --color-scheme classic -o $(logo_folder)/d$(i).png`);

        # do it for the reverse complement as well
        io = open(logo_folder*"/d$(i)_c.transfac", "w")
        println(io, "ID\t")
        println(io, "XX\t")
        println(io, "BF\t")
        println(io, "XX\t")
        println(io, "P0\tA\tC\tG\tT")
        q = Int.(floor.(reverse(ms.pfms[ind]) .* 1000));
        for j = 1:size(ms.pfms[ind],2)
            cur_rows = j < 10 ? string(Int(floor(j/10)))*"$j" : string(j);
            println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
        end
        println(io, "XX\t")
        close(io)            
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac -n 40 --errorbars NO -F png -X NO -Y NO --fineprint " " --resolution 600 -s large --fontsize 32 --color-scheme classic -o $(logo_folder)/d$(i)_c.png`);
    end
end

function save_pfms_as_transfac_noXY_weblogo(logo_folder::String, ms, sort_perm::Vector{Int})
    for i in sort_perm        
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i).transfac -n 40 --errorbars NO -F png -X NO -Y NO --fineprint " " --resolution 600 -s large --fontsize 32 --color-scheme classic -o $(logo_folder)/d$(i)_noxy.png`);      
        Base.run(`weblogo -D transfac -f $(logo_folder)/d$(i)_c.transfac -n 40 --errorbars NO -F png -X NO  -Y NO --fineprint " " --resolution 600 -s large --fontsize 32 --color-scheme classic -o $(logo_folder)/d$(i)_c_noxy.png`);
    end
end

# function get_html_template()



get_folder_names(target_folder::String) = 
    target_folder*"/"*logo_folder_name, target_folder*"/"*pics_folder_name

function make_folder_paths(target_folder, logo_folder, pics_folder)
    !isdir(target_folder) && (mkpath(target_folder);)
    !isdir(logo_folder) && (mkpath(logo_folder);)
    !isdir(pics_folder) && (mkpath(pics_folder);)
end

function render_main_summary_page(labels, 
                                  pvalues, 
                                  logos,
                                  target_folder,
                                  valid_alphas,
                                  activate_counts::Vector{Int},
                                  num_seqs
                                  )
    df = DataFrame(label=labels, eval=pvalues, logo=logos, counts=activate_counts);
    if length(valid_alphas) > 0
        out = Mustache.render(html_template_has_alpha_valid, 
                          target_folder=target_folder, 
                          valid_alphas="$valid_alphas",
                          num_alphas="$(length(valid_alphas)-1)",
                          min_alpha="$(valid_alphas[1])",
                          num_seq=num_seqs,
                          logo_folder=logo_folder_name, DF=df);
    else
        out = Mustache.render(html_template_no_alpha_valid, 
                          target_folder=target_folder, num_seq=num_seqs,
                          logo_folder=logo_folder_name, DF=df);
    end
    io = open(target_folder*"/summary.html", "w")
    print(io, out);
    close(io)
end

function save_result(ms, data,
                     target_folder, pval_cut_off;
                     alphas=[-1, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
                     )

    logo_folder, pics_folder = get_folder_names(target_folder);
    make_folder_paths(target_folder, logo_folder, pics_folder);
    p_vec, use_vec = 
        get_significant_motif_positions_pvec(ms, data, pval_cut_off);
    sort_perm = sortperm(p_vec); # sort it according to pvalues (small to big)
    p_vec_sort_perm, use_vec_sort_perm = p_vec[sort_perm], use_vec[sort_perm]
    pvalues = get_rounded_pval.(p_vec_sort_perm, use_vec_sort_perm);
    labels  = ["D$j" for j = 1:ms.num_motifs];
    logos   = ["d$(j)" for j = 1:ms.num_motifs];

    save_pfms_as_transfac(logo_folder, ms, sort_perm);
    save_pfms_as_transfac_noXY_weblogo(logo_folder, ms, sort_perm);

    # plot co-occurrences
    pwms_img   = [CairoMakie.load(target_folder*"/logos/d$(i)_noxy.png") for i = 1:ms.num_motifs];
    pwms_img_c = [CairoMakie.load(target_folder*"/logos/d$(i)_c_noxy.png") for i = 1:ms.num_motifs];
    valid_alphas, activate_counts = plot_cooccurrence(ms, 
                      data, 
                      sort_perm, 
                      use_vec, 
                      pwms_img, 
                      pwms_img_c, 
                      logo_folder; 
                      alphas=alphas
                      )

    render_main_summary_page(labels, 
                             pvalues, 
                             logos,
                             target_folder,
                             valid_alphas,
                             activate_counts,
                             data.N+data.N_test
                             ) 
    return pvalues, sort_perm
end