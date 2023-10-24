using Plots

L = 64
T = 64 #number of timesteps
p = 0.1

glider = [0 0 0 0 0 # 5x5 array containing a glider
          0 0 1 0 0
          0 0 0 1 0
          0 1 1 1 0
          0 0 0 0 0]

# gives moore neighbourhood of i,j as vector of tuples
neighbours(i::Integer, j::Integer) = [(i+1, j-1), (i+1, j), (i+1, j+1),
                                      (i,   j-1),           (i,   j+1),
                                      (i-1, j-1), (i-1, j), (i-1, j+1)]

periodic(i) = (i-1+L) % L + 1

function psi(M::Matrix, i, j)
    n = sum( [M[periodic(a), periodic(b)] for (a,b) in neighbours(i, j)] ) 
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

step(M::Matrix) = [psi(M, i, j) for i=1:L, j=1:L]

function main()
    color = cgrad([:white, :gray])
    M = zeros(Int8, L, L)
    #M = [Int(rand() < p) for i=1:L, j=1:L]
    M[1:5, 1:5] = glider
    mygif = @animate for t in 1:T
        heatmap(M, c = color)
        M = step(M)
    end

    gif(mygif, "Week06/game_of_life.gif", fps = 10)
end

main()