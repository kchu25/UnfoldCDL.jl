# checked for labeled data
function overlap_w_motif(data, pos, k)
    ol = 0
    for p in pos[k]
        ol += Cdlunroll.num_overlap(data.raw_data[p[2]].motif_where, p[1])
    end
    ol / (length(pos[k])*get_pos_len(pos[k]))
end

function olap_ratio(q1, q2)
    q1_a = get_pos_len(q1)*length(q1)
    q2_a = get_pos_len(q2)*length(q2)
    ol = 0
    for i in q1
        for j in q2
            if i[2] == j[2]
                ol += Cdlunroll.num_overlap(i[1],j[1])
                break
            end
        end
    end
    ol/(q1_a + q2_a - ol)
end

function chk_olap(new_pos)
    len_pos = length(new_pos)
    omat = zeros(Float64, (len_pos, len_pos))        
    for (ind1,i) in enumerate(new_pos)
        for ind2 = ind1+1:length(new_pos)
            if ind1 != ind2
                omat[ind1, ind2] = olap_ratio(i, new_pos[ind2])
            end
        end
    end
    return sort(omat[:], rev=true)
end
