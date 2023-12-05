using Plots
using Printf
using LinearAlgebra
using Dates

const dim = 3  # It is assumed that all simulations are in 3D.

###################################
# Structures
###################################

"Planet phase space at time t_n, and at some saved timesteps."
mutable struct Planet
    x::Vector{Float64}  # Position at time t_n.
    v::Vector{Float64}  # Velocity at time t_n, or t_n-1/2 if you're using a staggered time.
    a::Vector{Float64}  # Acceleration at time t_n.
    m::Float64  # Mass of planet.
    
    xs::Matrix{Float64}  # History of saved positions (note that not every timestep will be saved).
    vs::Matrix{Float64}  # History of saved velocities (note that not every timestep will be saved).

    Planet(x0::Vector{Float64}, v0::Vector{Float64}, m::Float64) = new(x0, v0, zeros(Float64, dim), m,
                                                                       Matrix{Float64}(undef,1,1), Matrix{Float64}(undef,1,1))
end

"Universe at time nΔt"
struct Universe
    planets::Vector{Planet}  # All planets (or other celestial bodies) in the universe.
    Δt::Float64  # Timestep.
    T::Float64  # Final time.
    ntot::Int64  # Number of timesteps.
    ns::Vector{Int64}  # Timesteps where position and velocity are saved.
    
    function Universe(planets::Vector{Planet}, Δt::Float64, T::Float64, nsave::Int64)
        # Get number of timesteps and correct final time.
        ntot = floor(Int64, T / Δt)
        T = ntot * Δt
        @printf "Final time of universe = %.2E\n" T        
        @printf "Number of steps = %.2E\n" ntot

        # Get timesteps to be saved.
        if nsave > ntot
            println("Warning: nsave > number of timesteps, so I will save every timestep, i.e. nsave = ntot")
            nsave = ntot
        end
        ns = round.(Int64, LinRange(1, ntot, nsave))
        
        # Allocate space for history of each planet.
        for planet in planets
            planet.xs = zeros(Float64, dim, nsave)
            planet.vs = zeros(Float64, dim, nsave)
        end

        # Create universe.
        new(planets, Δt, T, ntot, ns)
    end
end


###################################
# Functions
###################################

function simulate!(uni::Universe, integrator::Symbol; verbose = false)
    
    init_time = time()  # Measure time it takes to run simulation.
    N = length(uni.planets)  # N bodies of the N-body probem..

    # If leapfrog (with the staggered velocity), we need to shift v_0 to v_{-1/2}.
    if integrator == :leapfrog_staggered
        calculate_forces!(uni.planets)
        for planet in uni.planets
            planet.v .-= planet.a * uni.Δt/2
        end
    end

    # Run the simulation.
    idx_save = 1
    for n = 1:uni.ntot
        step!(uni, integrator)
        # Save current timestep.
        if n in uni.ns
            if verbose
                @printf "It took %.2F seconds to run %.2F %% of the simulation.\n" time()-init_time n/uni.ntot*100
            end
            for planet in uni.planets
                planet.xs[:,idx_save] .= planet.x
                planet.vs[:,idx_save] .= planet.v
            end
            idx_save += 1
        end
    end
end


function step!(uni::Universe, integrator::Symbol)
    Δt = uni.Δt
    calculate_forces!(uni.planets)

    if integrator == :euler
        for planet in uni.planets
            planet.x .+= planet.v * Δt  # x_{n+1} = x_n + v_n * Δt .
            planet.v .+= planet.a * Δt  # v_{n+1} = v_n + a_n * Δt .
        end
        
    elseif integrator == :leapfrog_staggered
        for planet in uni.planets
            planet.v .+= planet.a * Δt  # v_{n+1/2} = v_{n-1/2} + a_n * Δt .
            planet.x .+= planet.v * Δt  # x_{n+1} = x_n + v_{n+1/2} * Δt .
        end

    elseif integrator == :leapfrog_integer
        for planet in uni.planets
            planet.v .+= planet.a * Δt/2  # v_{n+1/2} = v_n + a_n * Δt/2 .
            planet.x .+= planet.v * Δt  # x_{n+1} = x_n + v_{n+1/2} * Δt .
        end
        calculate_forces!(uni.planets)
        for planet in uni.planets
            planet.v .+= planet.a * Δt/2  # v_{n+1} = v_{n+1/2} + a_{n+1} * Δt/2 .
        end

    elseif integrator == :RK4
        # Get current state.
        ms = [planet.m for planet in uni.planets]
        state = [[planet.x, planet.v] for planet in uni.planets]
        # Get intermediate increments.
        k1 = [[planet.v, planet.a] for planet in uni.planets]
        k2 = calculate_forces(state .+ Δt/2 * k1, ms)
        k3 = calculate_forces(state .+ Δt/2 * k2, ms)
        k4 = calculate_forces(state .+ Δt * k3, ms)
        # Update state.
        state .+= Δt/6 * (k1 .+ 2*k2 .+ 2*k3 .+ k4)
        for (p,planet) in enumerate(uni.planets)
            planet.x .= state[p][1]
            planet.v .= state[p][2]
        end

    else
        println("Unknown integrator method ", integrator)
    end

end


function calculate_forces!(planets::Vector{Planet})
    N = length(planets)

    # Reset accelerations to zero.
    for planet in planets
        planet.a .= zero(Float64)
    end

    # Now calculate forces.
    for i = 1:N-1
        for j = i+1:N
            f = get_force(planets[i].x, planets[j].x)
            planets[i].a .+= f * planets[j].m
            planets[j].a .-= f * planets[i].m
        end
    end
end

function calculate_forces(state::Vector{Vector{Vector{Float64}}}, ms::Vector{Float64})
    N = length(state)
    planets = [Planet(state[i][1], state[i][2], ms[i]) for i = 1:N]
    calculate_forces!(planets)
    return [[planet.v, planet.a] for planet in planets]
end

function get_force(x1::Vector{Float64}, x2::Vector{Float64})
    # The formula for the force is e*G*m1*m2/(x1-x2)^2 where e is a unit vector in the direction of x1-x2.
    # But since the acceleration is force/mass, we can just leave the masses out of the equation and multiply them later.
    # Thus we can just calculate f = e*G/(x1-x2)^2, or f = G*d/d^3 with d = x1-x2.
    d = x1 .- x2
    return -d .* G / norm(d)^3
end


function circle_points(r::Float64; x::Float64 = 0.0, y::Float64 =0.0)
    θ = LinRange(0,2π,100)
    return x .+ r * cos.(θ), y .+ r * sin.(θ)
end


ekin(ms::Vector{Float64}, vs::Vector{Vector{Float64}}) = sum([0.5 * ms[p] * norm(vs[p])^2 for p = 1:length(ms)])
ekin(planets::Vector{Planet}, n::Int64) = ekin([planet.m for planet in planets], [planet.vs[:,n] for planet in planets])


function epot(ms::Vector{Float64}, xs::Vector{Vector{Float64}})
    N = length(ms)
    E = 0
    for i = 1:N-1
        for j = i+1:N
            d = norm(xs[i] .- xs[j])
            E += G * ms[i] * ms[j] / d
        end
    end
    return E
end
epot(planets::Vector{Planet}, n::Int64) = epot([planet.m for planet in planets], [planet.xs[:,n] for planet in planets])

etot(planets::Vector{Planet}, n::Int64) = ekin(planets, n) + epot(planets, n)