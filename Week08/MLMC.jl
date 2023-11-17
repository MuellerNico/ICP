using Plots
using Distributions

ϵ = 10e-1
σ = 1
p = 1
ξ = Normal(0, 1)

function sde(x0, l, T, N) 
    x = x0
    dt = T / 2^l # 2ˡ timesteps
    for i=1:2^l
        x = x - x^p*dt + sqrt(2*dt)*σ*rand(ξ)
    end
    return x #x(t=T)
end

function mlmc()

end