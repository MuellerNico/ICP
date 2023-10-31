using Plots

f(x) = exp(-x^2)
C = 1.0 / (-exp(-1) + 1) # constant to make integral of g(x) 1
g(x) = exp(-x) * C
U(x) = 1

# generates point distribution based on g with the rejection method
function get_distribution(g::Function, N::Integer, RNG=rand)::Vector{Float64}
    dist = zeros(N)
    for i=1:N
        while true
            x, y = RNG(2)
            if y < g(x) # accept
                dist[i] = x
                break
            end
        end
    end
    return dist
end

# MC integration of f using the point distribution dist
function integrateMC(f::Function, dist::Vector{Float64})::Float64
    I = 0.0 # Integral
    for x in dist
        I += f(x)
    end
    return I / length(dist)
end

# quadrature for "exact" solution using high number of equidistant nodes (N ~ 100'000'000)
function trapezoidal(f::Function, N::Integer)::Float64
    h = 1.0 / (N-1) # step size
    I = 0.5 * f(0) # first node
    for x in range(h, 1-h, N-2) 
        I += f(x) # inner nodes
    end
    I += 0.5 * f(1) # last node
    I *= h
    return I
end

function plot_convergence()
    m = 8 # number of measurements. 8 is best
    N = [10^i for i=0:m-1]
    println("calculating exact value...")
    exact = trapezoidal(f, 100000000) # 10^8 nodes
    errUniform = zeros(m)
    errImportance = zeros(m)

    for i=1:m
        println("n = $(N[i])...")
        errUniform[i] = abs(integrateMC(f, get_distribution(U, N[i])) - exact)
        errImportance[i] = abs(integrateMC(x -> f(x)/g(x), get_distribution(g, N[i])) - exact)
    end

    println("creating plot...")
    plot(N, errUniform, label="uniform sampling")
    plot!(N, errImportance, label="importance sampling")
    plot!(xaxis=:log, yaxis=:log, xlabel="N (number of samples)", ylabel="absolute error")
    plot!(title="Error convergence of Monte Carlo integration")
    png("Week07/error_convergence")
    println("done.")
end

plot_convergence()

histogram(get_distribution(g, 10000), title="sampling using g(x)")
plot!([0:0.05:1], g, ylabel="g(x)")
png("Week07/g_dist")