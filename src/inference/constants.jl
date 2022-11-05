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

const vec_tup_t = NamedTuple{(:at, :f, :len, :f1_prev, :f2_prev, :f1_prev_loc, :f2_prev_loc), 
                      NTuple{7, Int}}

const pos_t = Tuple{UnitRange{Int}, Int};

const topo_t = NamedTuple{(:f1, :f2, :f3, :d12, :d13), NTuple{5, Int}};

const ndict_c_t = Dict{Tuple{Int,Int,Int,Int}, Bool}; # ndict's nested (child) dict type

const ndict_t = Dict{Int, ndict_c_t};



dna_comp(a::String) = join([DNA_complement[l] for l in reverse(a)]);