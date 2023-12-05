using Plots
using Distributions

L = 20 # domain length
T = 50 # timespan
num_timesteps = 1000 # num timesteps

lambda = 1/2 #dt/dx

dt = T / num_timesteps # timestep size
dx = dt / lambda
num_gridpoints = Int(L / dx) # num grid points

@show dt / dx
@show num_gridpoints
@show num_timesteps

t0 = 3
sigma = 1

E = zeros(num_gridpoints, num_timesteps)
B = zeros(num_gridpoints, num_timesteps)
N = Normal(L/2,5)

# BC
B[1,:] .= 0.0
B[end,:] .= 0.0
# IC
#E[:,1] = rand(N, num_gridpoints)

function integrate()
    for n = 1:num_timesteps-1
        for i = 1:num_gridpoints-1
            E[i,n+1] = E[i,n] - lambda*(B[i+1,n] - B[i,n])
        end
        # source
        E[Int(num_gridpoints/2),n+1] = exp(-((n*dt-t0)/sigma)^2)

        for i = 1:num_gridpoints-2
            B[i+1,n+1] = B[i+1,n] - lambda*(E[i+1,n+1] - E[i,n+1])
        end  
    end
end

function generate_gif()
    E_range = LinRange(1, num_gridpoints, num_gridpoints)
    B_range = LinRange(1, num_gridpoints, num_gridpoints) .+ 0.5
    anim = @animate for n = 1:num_timesteps
        plot(xlabel="L", ylabel="amplitude", ylims=(-2,2))
        plot!(E_range, E[:,n], label="E")
        plot!(B_range, B[:,n], label="B")
    end
    gif(anim, "Week12/FDTD.gif", fps=15)
end

function field_heatmap()
    heatmap(E, xlabel="time", ylabel="i")
    png("Week12/E_field")
    heatmap(B, xlabel="time", ylabel="i")
    png("Week12/B_field")
end

function plot_source()
    plot(E[Int(num_gridpoints/2),:], xlabel="timestep n", ylabel="value of E at L/2")
    png("Week12/source")
end

integrate()
field_heatmap()
generate_gif()
plot_source()