#= customized rules
See: https://discourse.julialang.org/t/forwarddiff-and-zygote-cannot-automatically-differentiate-ad-function-from-c-n-to-r-that-uses-fft/52440/3
=#

ifft_ZᵀZ(z, iFzz, Fzz) = iFzz * z;
Zygote.@adjoint ifft_ZᵀZ(z, iFzz, Fzz) = 
    ifft_ZᵀZ(z, iFzz, Fzz), c̄ -> (1 ./ length(c̄) .* (Fzz * c̄), nothing,  nothing);

ifft_ZᵀS(z, iFzs, Fzs) = iFzs * z;
Zygote.@adjoint ifft_ZᵀS(z, iFzs, Fzs) = 
    ifft_ZᵀS(z, iFzs, Fzs), c̄ -> (1 ./ length(c̄) .* (Fzs * c̄), nothing,  nothing)

irfft_ZᵀZ(z, irFzz, rFzz) = irFzz * z;
Zygote.@adjoint irfft_ZᵀZ(z, irFzz, rFzz) = 
    irfft_ZᵀZ(z, irFzz, rFzz), c̄ -> (1 ./ length(c̄) .* (rFzz * c̄), nothing,  nothing);

irfft_ZᵀS(z, irFzs, rFzs) = irFzs * z;
Zygote.@adjoint irfft_ZᵀS(z, irFzs, rFzs) = 
    irfft_ZᵀS(z, irFzs, rFzs), c̄ -> (1 ./ length(c̄) .* (rFzs * c̄), nothing,  nothing)
