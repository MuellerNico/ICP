using LinearAlgebra
using Plots
using Distributions

# Domain = [0,1]
n = 20 # space discretization
dx = 1.0 / n
dt = 0.01
u = 1.0 # velocity field
sigma = 0.1
T = 2
fps = 10

initial() = [1/sqrt(2*pi*sigma)*exp(-(j*dx)^2/(2*sigma^2)) for j=1:n] # x = j*dx
initial_square() = [0.3 < j*dx < 0.7 ? 1 : 0 for j=1:n]


function FTBS()
    c = u*dt/dx
    M = zeros(n,n)
    M[diagind(M)] .= 1-c # diagonal entries
    for j=2:n  # set subdiagonal
        M[j,j-1] = c
    end
    M[1,n] = c #special entry top right. periodic BDC
    
    phi = initial_square()
    x = LinRange(0, 1, n)

    mygif = @animate for t = 0:dt:T
        phi = M * phi
        plot(x, phi, xaxis="x", yaxis="ϕ(x,t)", ylims=(-1,1), title="FTBS, n=$n, speed=$(fps*dt), t=$t")
    end

    gif(mygif, "Week09/advection_FTBS.gif", fps = fps)
end

function CTCS()
    c = u*dt/dx
    M = zeros(2n, 2n)
    M[n+1:2n, 1:n] += I # lower left identity
    M[1:n, n+1:2n] += I # upper right identity
    for j=1:n
        M[j, (j)%n+1] = -c # superdiagonal
        M[j, (j-2+n)%n+1] = c #subdiagonal
    end
    heatmap(M, yflip = true)
    png("matrixheatmap")

    phi = zeros(2n)
    phi[1:n] = initial_square()
    phi[n+1:2n] = initial_square() # not correct. use FTCS for first timestep

    x = LinRange(0, 1, n)
    mygif = @animate for t = 0:dt:T
        phi = M * phi
        plot(x, phi[1:n], xaxis="x", yaxis="ϕ(x,t)", ylims=(-1,1), title="FTBS, n=$n, speed=$(fps*dt), t=$t")
        plot!(x, phi[n+1:2n])
    end
    gif(mygif, "Week09/advection_CTCS.gif", fps = fps)
end

CTCS()

