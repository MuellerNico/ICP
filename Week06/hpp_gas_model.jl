using Plots 

L = 128
T = 128
p = 0.02

Cell = Vector{Bool}    # data type
cell()::Cell = [0, 0, 0, 0]   # (north, east, south, west)
rand_cell(p)::Cell = [rand() < p for _=1:4]

# manipulates indices to apply periodic BDC
periodic(i) = (i-1+L) % L + 1 

neighbors(i, j) = [(periodic(x), periodic(y)) for (x,y) in [(i-1,j), (i,j+1), (i+1,j), (i,j-1)]] # [north, east, south, west]

display(M) = maximum.(M) # each cell with 1 or more particles is filled, rest empty
display_alt(M) = sum.(M) # alternative to show shades depending on density within cell

# I think I messed up the implementation of this part
function gaussian_dist()::Matrix{Cell}
    sigma = 5
    p(i,j) = exp(- ((i-L/2)^2 + (j-L/2)^2) / sigma^2)
    M = [cell() for i=1:L, j=1:L]
    for i=1:L, j=1:L
        for d=1:4
            if rand() < p(i,j)
                M[i,j][d] = 1
            end
        end
    end
    return M
end

# compute cell state after collision
function collision(c::Cell)::Cell
    if c == [1,0,1,0]
        [0,1,0,1]
    elseif c == [0,1,0,1]
        [1,0,1,0]
    else 
        c
    end
end

# compute cell state based on neighbours at t-1
function propagate_cell(M::Matrix{Cell}, i, j)::Cell
    c = cell()
    for d=1:4 # looop through directions 1-4 (north-west)
        (x,y) = neighbors(i, j)[d]
        c[d] = M[x,y][d]
    end
    return c
end

function main()
    color = cgrad([:white, :gray])
    #M = [rand_cell(p) for i=1:L, j=1:L]
    #M[16,1] = [0,0,0,1]
    #M[16,31] = [0,1,0,0] 
    M = gaussian_dist()

    mygif = @animate for t=1:T
        heatmap(display_alt(M), clims=(0,4), c=color, yflip=true)
        title!("HPP Gas Model. $T timesteps")
        M = collision.(M)
        M = [propagate_cell(M, i, j) for i=1:L, j=1:L]
    end
    gif(mygif, "Week06/hpp.gif", fps = 10)
end

main()

#= -- notes --
particles can "jump" over each other with this way of timestepping
=#