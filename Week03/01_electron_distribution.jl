using Roots
using Plots
using BenchmarkTools

# println(Threads.nthreads())

N = 10000 #use multiple of 1'000 
g(x) = 1/(4*pi) * (1 + 0.5*cos(x/2)) # normalized g(x)
I(x) = 1/(4*pi) * (x + sin(x/2)) # normalized I(x)

function generateRN(n)
    rn = zeros(n)
    for i in 1:n
        u = rand()
        f(x) = I(x) - u
        rn[i] = find_zero(f, (0,4*pi))
    end
    return rn
end

function generateRN_threading(n)
    numThreads = Threads.nthreads()
    #println("Running RN generation with $numThreads Threads...")
    rn = zeros(n)
    chunkSize = trunc(Int, n/numThreads)
    Threads.@threads for i in 1:numThreads
        for j in 1:chunkSize
            u = rand()
            f(x) = I(x) - u
            rn[(i-1)*chunkSize + j] = find_zero(f, (0,4*pi))
        end
    end
    return rn
end

function plot_rn()
    x = generateRN_threading(N)
    b_range = range(0, 4*pi, length=12)
    histogram(x, bins=b_range, normalize=:pdf, label="transformation method")
    plot!(g, 0, 4*pi, label="g(x)")
    xlabel!("x")
    ylabel!("P(x)")
    title!("N = $N")
    png("Week03/1D")
end

function benchmark_rn()
    pow = 5 # power of ten for max number
    s = zeros(pow) #speedup
    n = zeros(Int32, pow)
    for i in 1:pow
        n[i] = 10^i
        t1_test = generateRN(n[i])
        t1 = mean(t1_test)
        tp_test = generateRN_threading(n[i])
        tp = mean(tp_test)
        s[i] =  t1 / tp
    end
    println(n)
    scatter(n, s)
    plot!(xscale = :log10)
    xlabel!("n (number of particles)")
    ylabel!("S = T₁/Tₚ")
    title!("Speedup on $(Threads.nthreads()) cores")
    png("Week03/speedup")
end

benchmark_rn()