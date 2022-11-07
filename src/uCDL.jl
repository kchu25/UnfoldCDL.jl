module uCDL

using Flux, Zygote, CUDA, FFTW, LinearAlgebra, Random, 
      StatsBase, Distributions, StaticArrays, DoubleFloats,
      DataStructures, HypothesisTests, SimDNA, Mustache,
      DataFrames, Makie, CairoMakie, ColorSchemes, FastaLoader

export find_motif, find_motif_fasta_folder

include("constants.jl")
include("custom_Zygote_rules.jl")
include("utils.jl")
include("model_basic.jl")
include("opt.jl")
include("training.jl")
include("inference/constants.jl")
include("inference/triplet_esd_test.jl")
include("inference/triplet.jl")
include("inference/utils.jl")
include("inference/merging.jl")
include("inference/expansion.jl")
include("inference/inference.jl")
include("greedy_alignment/constants.jl")
include("greedy_alignment/PWM_Touzet.jl")
include("greedy_alignment/utils.jl")
include("greedy_alignment/expansion.jl")
include("greedy_alignment/scan.jl")
include("greedy_alignment/scan_gpu.jl")
include("greedy_alignment/fisher.jl")
include("greedy_alignment/greedy_alignment.jl")
include("performance_eval/html_template.jl")
include("performance_eval/save_cooccurrence.jl")
include("performance_eval/eval.jl")
include("performance_eval/eval_test_set.jl")
include("performance_eval/save.jl")
include("performance_eval/jaspar.jl")
include("find_motif.jl")

end
