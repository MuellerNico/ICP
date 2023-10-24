
using Plots

function new_lattice(n::Integer, p::Real)::Matrix{Int8}
    lat = zeros(Int8, n+2,n+2)
    lat[2:end-1, 2:end-1] = [Int8(rand() < p) for i=1:n, j=1:n] 
    return lat
end

function update!(G::Matrix{Int8})::Bool
    L = size(G)[1]-1
    stop::Bool = true
    for i=2:L, j=2:L    # loop over entire grid
        if G[i,j] == 2  # if burning, burn neighbours
            stop = false
            for (a,b) in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]
                if G[a,b] == 1
                    G[a,b] = 2
                end
            end
            G[i,j] = 3 #b urnt down
        end
    end
    return stop # true if no new trees got lit, fire stopped
end

function draw_forest(G::Matrix{Int8}, suffix::String="")
    color = cgrad([:white,:green,:red])
    heatmap(G, c=color, clims=(0, 3), yflip = true)
    png("Week05/forest" * suffix)
end

function percolated(G::Matrix{Int8})
    L = size(G)[1]-1
    for i=2:L
        if G[i,L] == 3
            return true
        end    
    end
    return false
end

# sets all burned trees (3) to 1 and everything else to 0
function clean!(G::Matrix{Int8})
    for i=1:size(G)[1], j=1:size(G)[2]
        G[i,j] = (G[i,j] == 3 ? 1 : 0)
    end
    return G
end

function compute_cluster(G::Matrix{Int8})
    L = size(G)[1]-1    
    for i=2:L 
        Gi = copy(G) 
        if Gi[i,2] == 1
            Gi[i,2] = 2 # try with this site
            while !update!(Gi)
            end
            if percolated(Gi)
                clean!(Gi)
                return (true, Gi, i)
            end
        end
    end
    return (false, G, 0)
end

M = new_lattice(32, 0.595)
#print(M)
draw_forest(M, "_before")
res = compute_cluster(M)
draw_forest(res[2], "_after")
if res[1] == false
    println("did not percolate")
else 
    println("percolated from idx $(res[3])")
end

