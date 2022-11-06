Base.@kwdef mutable struct HyperParam
    batch_size::Int = 300       # batch size for NN
    filter_len::Int = 6         # length of each filter
    f_len::Int = 6*4            # filter length for dummy encoded DNA sequences
    M::Int = 85                 # number of filters
    K_c::Int = 4                # iteration number for CSC
    K_d::Int = 1                # iteration number for CDL
    gamma_1::float_type=45.5    # regularization param for incoherence
    pool_mask_len::Int=40
    pool_stride::Int=4
end

function randomly_initialize_filters(; 
                       dim=4,
                       rng=Random.GLOBAL_RNG, 
                       repeats=5, 
                       how_many_filters=10,
                       float_type=Float16)    
    arr = zeros(float_type,(dim+1,
                            repeats,
                            1,
                            how_many_filters));    
    @inbounds for i = 1:repeats, j = 1:how_many_filters
        unif = rand(rng, dim-1);
        arr[2:dim,i,1,j] .= sort(unif);
    end
    @inbounds arr[dim+1,:,:,:] .= 1;
    return reshape(diff(arr, dims=1), (dim*repeats,1,how_many_filters));
end

struct CDL
    D1::CuArray{float_type, 3}
    D2::CuArray{float_type, 3}
    Dfinal::CuArray{float_type, 3}
    eta::Array{float_type, 1}
    lambda::Array{float_type, 1}
    mu::Array{float_type, 1}

    function CDL(hp, K_d, K_c_p1; η₁=0.0005, η₂=-0.001)    
        # K_d: number of iterations for CDL
        # K_c_p1: number of iterations for CSC plus 1 (for the initial iteration)
        # η₁,η₂: factor to make the number small
        D1 = cu(randomly_initialize_filters(
                                   repeats=hp.filter_len, 
                                   how_many_filters=hp.M,
                                   float_type=float_type));    
        D2 = cu(randomly_initialize_filters(
                        repeats=hp.filter_len, 
                        how_many_filters=hp.M,
                        float_type=float_type));    
        Dfinal = cu(randomly_initialize_filters(
                        repeats=hp.filter_len, 
                        how_many_filters=hp.M,
                        float_type=float_type));    
        D1 = sqrt.(D1); D2 = sqrt.(D2); D1 = sqrt.(Dfinal); 
        eta    = float_type(η₁)*rand(float_type, K_c_p1);
        lambda = float_type(η₁)*rand(float_type, K_c_p1);
        mu     = float_type(η₂)*rand(float_type, K_d);
        new(D1, D2, Dfinal, eta, lambda, mu)
    end
end

Flux.@functor CDL
