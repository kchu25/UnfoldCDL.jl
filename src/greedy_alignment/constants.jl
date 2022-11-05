############# for PWM_Touzet.jl ################
const _granularity_ = 1e-1; # initial granularity for score2pvalue and pval2score
const _k_ = 100; # decreasing factor for finer granularity in each iteration
const _bg_ = [.25,.25,.25,.25]; # default background
################################################

const pvalue_pwm_thresh = 0.00027
const pvalue_alpha_fisher = 1e-5
const atcg_comp = Dict{Char, Char}('a'=>'t', 'A'=>'T', 
                                   'c'=>'g', 'C'=>'G',
                                   'g'=>'c', 'G'=>'C',
                                   't'=>'a', 'T'=>'A')
# kernel thread-block setup
const threads_1d = 512;
const threads_2d = 32;
const threads_3d = 10;
const ker_1d = threads_1d;
const ker_2d = (threads_2d, threads_2d);
const ker_3d = (threads_3d, threads_3d, threads_3d);

b_size_1d(X) = ceil.(Int, size(X) ./ threads_1d)
b_size_2d(X) = ceil.(Int, size(X) ./ threads_2d)
b_size_3d(X) = ceil.(Int, size(X) ./ threads_3d)