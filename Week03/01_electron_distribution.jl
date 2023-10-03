using Roots
using Plots

# println(Threads.nthreads())

N = 10000
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
    rn = zeros(n)
    Threads.@threads for i in 1:n
        u = rand()
        f(x) = I(x) - u
        rn[i] = find_zero(f, (0,4*pi))
    end
    return rn
end

function benchmark()
    for n in range(100, 10000, 500)
        
        #generateRN(n)
    end
end

#scatter(rn)
x = generateRN(N)
b_range = range(0, 4*pi, length=12)
histogram(x, bins=b_range, normalize=:pdf)
plot!(g, 0, 4*pi)
xlabel!("x")
ylabel!("P(x)")
png("Week03/1D")
