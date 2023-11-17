using Plots

L = 100
R = 5

dist((x1,y1,z1), (x2,y2,z2)) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

# tries to place a random sphere. returns if successful or not
function place!(spheres::Vector)::Bool  
    p = Tuple(R .+ (L-2*R) * rand(3)) # random 3D coords between R and L-R in both x and y 
    for q in spheres
        if dist(p, q) < 2*R
            return false
        end
    end
    push!(spheres, p)
    return true
end

# creates a new configuration of n particles
function config(n::Integer)::Vector{Tuple}
    spheres = []
    limit = 200 # max number of placement attempts per particle
    for i=1:n
        attempts = 0
        while !place!(spheres)
            attempts += 1
            if attempts == limit
                println("ERROR: Could not place spheres")
                return nothing
            end
        end
    end
    return spheres
end

# ̅dₖ
function dk(spheres::Vector{Tuple})::Float64
    n = length(spheres)
    res = 0.0
    for i=1:n
        for j=i+1:n
            res += dist(spheres[i], spheres[j])
        end
    end
    return 2 / (n*(n-1)) * res
end

avg_dist2(M, n) = sum([dk(config(n)) for _=1:M]) / M

# average distance over M configuration with n particles each
function avg_dist(M, n)::Float64
    res = 0.0
    for i=1:M
        spheres = config(n)
        res += dk(spheres)
    end
    return res / M
end

function main()
    M = 1:10:250
    n = 5
    d = [avg_dist(m, n) for m in M]
    plot(M, d, xlabel="M (number of configurations)", ylabel="<d> (average distance)", ylims=(0, L))
    plot!(title="Integral convergence. L=$L, R=$R, n=$n")
    png("Week07/m-convergence")

    m = 1
    N = 2:10:300
    d = [avg_dist(m, n) for n in N]
    plot(N, d, xlabel="N (number of particles)", ylabel="<d> (average distance)", ylims=(0,L))
    plot!(title="n convergence. L=$L, R=$R, M=$m")
    png("Week07/n-convergence")
end

main()
