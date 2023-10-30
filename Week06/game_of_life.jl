using Plots

L = 64
T = 256 #number of timesteps
p = 0.15

glider = [0 0 0 0 0 # 5x5 array containing a glider
          0 0 1 0 0
          0 0 0 1 0
          0 1 1 1 0
          0 0 0 0 0]

# gives moore neighborhood of i,j as vector of tuples
neighbors(i::Integer, j::Integer) = [(i+1, j-1), (i+1, j), (i+1, j+1),
                                     (i,   j-1),           (i,   j+1),
                                     (i-1, j-1), (i-1, j), (i-1, j+1)]

# manipulates indices to apply periodic BDC
periodic(i) = (i-1+L) % L + 1 

step(M::Matrix) = [psi(M, i, j) for i=1:L, j=1:L]

function psi(M::Matrix, i, j)
    n = sum( [M[periodic(a), periodic(b)] for (a,b) in neighbors(i, j)] ) 
    if n == 0 || n == 1 #isolation
        0
    elseif n == 2 
        M[i,j]
    elseif n == 3 #birth
        1
    elseif n >= 4 #overpopulation
        0
    end
end

function insert_glider!(M::Matrix, i, j)
    M[i:i+4, j:j+4] = glider
end

function main()
    color = cgrad([:white, :gray])
    M = zeros(Int8, L, L)
    M = [Int(rand() < p) for i=1:L, j=1:L] # creating a random lattice to begin with
    #insert_glider!(M,1,1)
    #insert_glider!(M,20,20)
    #insert_glider!(M,30,4)
    mygif = @animate for t in 1:T
        heatmap(M, title = "Game of Life. $T timesteps", c = color)
        M = step(M) # timestep
    end

    gif(mygif, "Week06/game_of_life.gif", fps = 10)
end

main()