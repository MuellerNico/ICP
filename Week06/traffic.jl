using Plots

L = 32
T = 32
p = 0.7 # probability that initial square is occupied

@enum BDC begin
    dirichlet = 1
    periodic = 2
end

function simulate(bdc::BDC) 
    M = zeros(Int8, T, L+2) # each row i represents the road at timestep i
    v = [Int8(rand() < p) for _=1:L] #initial state
    M[1,2:L+1] = v # BDC: have extra zero row and column on both sides

    for t=2:T #timesteps
        for i = 2:L+1 #squares
            me = M[t-1, i]
            left = M[t-1, i-1]
            right = M[t-1, i+1]
            if bdc == periodic
                left =  M[t-1, (i-1-2+L) % L + 2] # i-1 (left) offset by -2 to get correct modulo
                right = M[t-1, (i+1-2)   % L + 2] # 1+1 (right) offset by -2
            end
            # rule 184
            if me == 0 && left == 1
                M[t,i] = 1
            elseif me == 1 && right == 1
                M[t,i] = 1
            else
                M[t,i] = 0
            end
        end
    end
    return M
end

function draw_traffic(M)
    color = cgrad([:white,:red])
    heatmap(M[:,2:end-1], c = color)
    ex, ey = range(0.5,L-0.5,L), range(0.5,T-0.5,T)
    vline!(ex, c=:black)
    hline!(ey, c=:black)
    plot!(xlabel="i", ylabel="t")
    plot!(title = "Traffic rule 184, L = $L")
    png("Week06/traffic")
end

draw_traffic(simulate(periodic))

#=
-- Conclusions --
Pattern resembles a checker board. With L/2 cars it's a perfect pattern. With less there are some gaps. 
With more cars congestion is created as some form of "back log" and creates a filled stripe along the spacetime plot. 

Phase transition at p=0.5
=#
