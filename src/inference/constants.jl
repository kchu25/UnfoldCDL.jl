const atcg2dummy = Dict{Char, Vector{Float32}}('a'=>[1,0,0,0], 'A'=>[1,0,0,0],
                                               'c'=>[0,1,0,0], 'C'=>[0,1,0,0],
                                               'g'=>[0,0,1,0], 'G'=>[0,0,1,0],
                                               't'=>[0,0,0,1], 'T'=>[0,0,0,1],
                                               'z'=>[0,0,0,0])

const vec_tup_t_ = Vector{Tuple{UnitRange{Int}, Int, Bool}};

const DNA_complement = Dict('A'=>'T', 'a'=>'t',
                      'C'=>'G', 'c'=>'g',
                      'G'=>'C', 'g'=>'c',
                      'T'=>'A', 't'=>'a');

dna_comp(a::String) = join([DNA_complement[l] for l in reverse(a)]);