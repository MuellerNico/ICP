using Distributed
using Plots
using Distributions

g = 9.81 # [m/s²] gravity
nu = 10e-6 # [m²/s] kinematic viscosity of water at 20°C
b = 0.47 # viscous drag coefficient (of a sphere)

L = 1.0 # domain length
N = 100 # number of grid points in both directions
T = 4.0 # duration
M = 400 # number of time steps
FPS = 30 # fps for the animation

@show dx = L/N
@show dt = T/M

r2(i,j) = (abs(i-N/2)*dx)^2 + (abs(j-N/2)*dx)^2 # calculate r² from the center
    
function initial_h()
    A = 20*dx # amplitude
    sigma = 0.2
    h = zeros(N+2,N+2) # one layer of ghost cells for boundary
    h[2:end-1,2:end-1] = [A * exp(- r2(i,j) / sigma^2) for i=1:N, j=1:N]
    surface(h, zlims=(0,1))
    png("Week10/initial_h")
    return h
end

# gradient of f in x direction
grad_x(f, i, j) = (f[i+1,j] - f[i-1,j]) / (2*dx)

# gradient of f in y direction
grad_y(f, i, j) = (f[i,j+1] - f[i,j-1]) / (2*dx)

# 2nd derivative in x direction
second_derivative_x(f, i, j) = ( f[i-1,j] - 2*f[i,j] + f[i+1,j] ) / (dx^2)

# 2nd derivative in y direction
second_derivative_y(f, i, j) = ( f[i,j-1] - 2*f[i,j] + f[i,j+1] ) / (dx^2)

function apply_PBDC!(M, m)
    M[1, 2:end-1]   = M[end-1, 2:end-1]
    M[end, 2:end-1] = M[2, 2:end-1]
    M[2:end-1, 1]   = M[2:end-1, end-1]
    M[2:end-1, end] = M[2:end-1, 2]
end


function main()

    velx = zeros(M, N+2, N+2)
    vely = zeros(M, N+2, N+2)
    sol = zeros(M, N+2, N+2)
    sol[1,:,:] = initial_h()

    @show size(velx)
    @show size(vely)
    @show size(sol)

    @show r2(50,50)
    @show r2(51,50)

    for m in 1:M-1
        t = m*dt
        if m%10 == 0
            println("m = $m / $M")
        end
        for i=2:N+1, j=2:N+1
            #@show (i,j) 
            #dh = - gradx(h.*vx,i,j) - grady(h.*vy,i,j)
            #dvx = -g*gradx(h,i,j) - (vx[i,j] * gradx(vx,i,j) + vy[i,j] * grady(vx,i,j)) - b*vx[i,j] + nu*(laplacex(vx,i,j) + laplacey(vx,i,j))
            #dvy = -g*grady(h,i,j) - (vx[i,j] * gradx(vy,i,j) + vy[i,j] * grady(vy,i,j)) - b*vy[i,j] + nu*(laplacex(vy,i,j) + laplacey(vy,i,j))
            
            vx = @view velx[m,:,:]
            vy = @view vely[m,:,:]
            h = @view sol[m,:,:]
            
            #apply_PBDC!(sol)
            #apply_PBDC!(velx)
            #apply_PBDC!(vely)

            dh_dt = -(grad_x(h,i,j)*vx[i,j] + grad_x(vx,i,j)*h[i,j] + grad_y(h,i,j)*vy[i,j] + grad_y(vy,i,j)*h[i,j])

            gravity = -g * [grad_x(h,i,j), grad_y(h,i,j)]
            velocity_grad_x = -(vx[i,j]*grad_x(vx,i,j) + vy[i,j]*grad_y(vx,i,j))
            velocity_grad_y = -(vx[i,j]*grad_x(vy,i,j) + vy[i,j]*grad_y(vy,i,j))
            velocity_grad = [velocity_grad_x, velocity_grad_y]
            drag = -[b*vx[i,j], b*vy[i,j]]
            kinematic_visc = nu*[second_derivative_x(vx,i,j), second_derivative_y(vy,i,j)] # this is just component wise 1D second derivative

            dvx_dt::Float64 = gravity[1] + velocity_grad[1] + drag[1] + #kinematic_visc[1]
            dvy_dt::Float64 = gravity[2] + velocity_grad[2] + drag[2] + #kinematic_visc[2]

            sol[m+1,i,j]   = h[i,j] + dt * dh_dt
            velx[m+1,i,j]  = vx[i,j] + dt * dvx_dt
            vely[m+1,i,j]  = vy[i,j] + dt * dvy_dt
            
            if sol[m+1,i,j] < 0 || sol[m+1,i,j] == NaN
                sol[m+1,i,j] = 0
            elseif sol[m+1,i,j] > 1
                sol[m+1,i,j] = 1
            end
        end
        # BDC
        #plot(LinRange(0,1,N+2), h[Int(N/2),:])
    end
    println("finished simulation.")
    println("plotting surface...")

    max_height = maximum(sol[1,:,:])
    max_velocity = maximum(velx)

    water_surface = @animate for m=1:M
        plot(surface(sol[m,:,:]), zlims=(0,max_height), title="t=$(m*dt)")
    end

    gif(water_surface, "Week10/solution.gif", fps = FPS)
    println("plotting velocity field...")

    vel_field = @animate for m=1:M
        abs_vel = [sqrt(velx[m,i,j]^2 + vely[m,i,j]^2) for i=1:N+2, j=1:N+2]
        heatmap(abs_vel, clims=(0,max_velocity), title="t=$(m*dt)")
    end

    gif(vel_field, "Week10/velocity_field.gif", fps = FPS)
    println("done.")
end

main()
