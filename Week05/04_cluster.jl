using Plots

function new_lattice(n::Integer, p::Real)::Matrix{Int32}
    lat = zeros(Int8, n+2,n+2)
    lat[2:end-1, 2:end-1] = [Int8(rand() < p) for i=1:n, j=1:n] 
    return lat
end

function recursive_lookup(M::Vector, k::Integer)::Integer
    #println("k = $k, M = $M") #debug
    if k == 0
        return 0
    elseif M[k] < 0
        return recursive_lookup(M, -M[k])
    else
        return k
    end
end

function hoshen_kopelman!(lat::Matrix)::Vector  #returns vector M of size L*L (ex 16*16=254)
    L1, L2 = size(lat)[1]-1, size(lat)[2]-1
    M = zeros(Int32, (L1-1)*(L2-1))  # index aka cluster size
    k = 2
    for j=2:L2
        for i=2:L1
            label = lat[i,j]
            left = recursive_lookup(M, lat[i,j-1])
            above = recursive_lookup(M, lat[i-1,j])
            if label == 1 # if populated
                if left == 0 && above == 0
                    lat[i,j] = k
                    M[k] = 1
                    k += 1
                elseif left != 0 && above == 0 # only left neighbour
                    lat[i,j] = left
                    M[left] += 1
                elseif left == 0 && above != 0 # only right neighbour
                    lat[i,j] = above
                    M[above] += 1
                elseif left == above
                    lat[i,j] = above
                    M[above] += 1
                else    # 2 neighbours => join
                    lat[i,j] = above
                    M[above] += 1 + M[left]
                    M[left] = -above
                end
            end
        end
    end
    for i=2:L1  # uniform coloring of clusters
        for j=2:L2
            lat[i,j] = recursive_lookup(M, lat[i,j])
        end
    end

    return M
end

# returns a vector where hist[x] returns the number of clusters of size x
function cluster_distribution(M::Vector)
    hist = zeros(length(M))
    for x in M
        if x > 0
            hist[x] += 1
        end
    end  
    # println("hist = $hist") #debug
    return hist  
end

function plot_cluster_size_dist()
    num_samples = 1000 #samples to average 
    p_steps = 10
    L = 16
    N = L * L
    Ns = zeros(N)
    np = []
    chip = []
    for p in range(0.1, 1, p_steps)
        for i=1:num_samples
            lat = new_lattice(L, p)
            M = hoshen_kopelman!(lat)
            Ns += cluster_distribution(M) # histogram
        end
        ns = 1/N * Ns
        push!(np, copy(ns))
        # calculate 2nd moment
        chi = 0
        idx_largest_cluster = maximum(findall(x -> x>0, Ns))
        N_clusters = sum(Ns)
        for s=1:idx_largest_cluster-1
            chi += s^2 * Ns[s] / N_clusters
        end
        push!(chip, chi)
    end
    plot(np, title="Cluster Distribution with $L x $L, $num_samples samples for all $p_steps discrete p valuesfor all $p_steps discrete p values")
    png("Week05/cluster-distribution")
    plot(chip, xlabel="p*10", ylabel="ᵪ", title="Chi (2nd moment) with $L x $L, $num_samples samples for all $p_steps discrete p values")
    png("Week05/second-moment")
end

function quick_test()
    L = 16
    p = 0.4
    lat = new_lattice(L, p)
    heatmap(lat, yflip = true)
    png("Week05/test-before")

    M = hoshen_kopelman!(lat)

    println("M ($(size(M))) = $M")
    heatmap(lat, yflip = true)
    png("Week05/test-after")
end

quick_test()
plot_cluster_size_dist()