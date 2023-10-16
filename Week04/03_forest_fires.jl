using Plots

L = 32
#=
function burn!(x::Int8)
    if x == 1
        x += 1
    end
end
=#
function update!(G::Matrix{Int8})::Bool
    stop::Bool = true
    for i=1:L, j=1:L    # loop over entire grid
        if G[i,j] == 2  # if burning, burn neighbours
            stop = false
            if i < L && G[i+1,j] == 1 #TODO move this to different funtion
                G[i+1,j] = 2
            end
            if i > 1 && G[i-1,j] == 1
                G[i-1,j] = 2
            end
            if j < L && G[i,j+1] == 1
                G[i,j+1] = 2
            end
            if j > 1 && G[i,j-1] == 1
                G[i,j-1] = 2
            end
            G[i,j] = 3 #burnt down
        end
    end
    return stop # true if no new trees got lit, fire stopped
end

function new_lattice(p)::Matrix{Int8}
    return [rand() < p ? 1 : 0 for i=1:L, j=1:L] 
end

function draw_forest(G::Matrix{Int8}, suffix::String="")
    color = cgrad([:white,:green,:red])
    heatmap(G, c=color, clims=(0, 2))
    png("Week04/forest" * suffix)
end

function percolating!(G::Matrix{Int8})
    # burn first column
    for i=1:L  
        if G[i,1] == 1
            G[i,1] = 2
        end
    end
    # update until finished
    while !update!(G)
    end
    # return true if any tree on other end is burned down
    for i=1:L
        if G[i,L] == 3
            return clustersize(G, i)
        end
    end
    return 0
end

function clustersize(G::Matrix{Int8}, idx)
    
end


M = new_lattice(0.69)
#print(M)
draw_forest(M, "_before")
res = percolating!(M)
println("Size of percolating cluster: $res")
draw_forest(M, "_after")

