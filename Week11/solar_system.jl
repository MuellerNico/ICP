using Plots
using LinearAlgebra
using Printf

# current issue: F_res probably applied wrongly. 

G = 6.674e-11 # [N*m²/kg²] Gravitational constant
JUPITER_MULTIPLIER = 1.0 # scale jupiters mass. with 10³ mars bounces from sun to jupiter, back, then gets lost.
M_SUN = 1988500E+24 # mass of the sun (used for E_pot calculation)

struct Planet
    name::String
    m::Float64 # [kg] Mass
    x::Vector{Vector{Float64}} # [m] num_timesteps x 2
    v::Vector{Vector{Float64}} # [m/s] num_timesteps x 2
    function Planet(name, m, x0::Vector{Float64}, v12::Vector{Float64}, num_timesteps::Integer)
        x = [zeros(2) for t=1:num_timesteps]
        v = [zeros(2) for t=1:num_timesteps]
        x[1] = x0 # initial condition
        v[1] = v12 # initial condition v_1/2
        return new(name, m, x, v)
    end
end

# TODO: move data to input file and read it here
function initial_config(num_timesteps::Integer)::Vector{Planet}
    earth = Planet("Earth", 5.97219E+24, 
                  [ 6.082184600127006E+07, 1.335381239845273E+08] * 10^3, # convert km to m
                  [-2.749517400867114E+01, 1.239113504931076E+01] * 10^3,
                  num_timesteps)

    sun = Planet("Sun", 1988500E+24,
                [1.216118823721544E+06, -3.980787627510761E+05] * 10^3,
                [0.0, 0.0], #[7.965896870326740E-03, -1.275819193112228E-02] * 10^3,
                num_timesteps)

    jupiter = Planet("Jupiter", 1898.18722E+24 * JUPITER_MULTIPLIER,
                    [ 5.484660608236390E+08, 5.019708824127624E+08] * 10^3,
                    [-8.965445298768357E+00, 1.025733727525911E+01] * 10^3,
                    num_timesteps)

    mars = Planet("Mars", 6.4171E+23,
                 [-1.141105292065805E+08, -1.986049892445633E+08] * 10^3,
                 [ 2.197326254440826E+01, -9.929422734484632E+00] * 10^3,
                 num_timesteps)

    neptune = Planet("Neptune", 102.409E+24,
                    [4.462373023571321E+09, -2.843814409809310E+08] * 10^3,
                    [3.099783650688389E-01,  5.456941175121609E+00] * 10^3,
                    num_timesteps)

    planets = [sun, earth, jupiter, mars]
    return planets
end

# returns force vector acting on planet 1
function newton_gravity(m1, m2, x1, x2)
    return -G * m1 * m2 / norm(x1 - x2)^2 * normalize(x1 - x2)
end

function resulting_force(p::Planet, planets::Vector{Planet}, t::Integer)
    return sum(newton_gravity(p.m, p2.m, p.x[t], p2.x[t]) for p2 in planets if p2 != p)
end

function integrate_leapfrog!(planets::Vector{Planet}, num_timesteps::Integer, dt)
    for t=1:num_timesteps-1

        if t%(num_timesteps/10) == 0
            println("t = $t / $num_timesteps")
        end
        
        for p in planets
            if p.name != "Sun" #don't update sun. seems to break otherwise
                p.x[t+1] = p.x[t] + p.v[t] * dt
                F_res = resulting_force(p, planets, t+1)
                acc = 1/p.m * F_res
                p.v[t+1] = p.v[t] + acc * dt
            end
        end
    end
end

function animate_orbit(planets, duration)
    FPS = 15
    num_timesteps = length(planets[1].x)
    num_frames = duration * FPS
    gif_step = Int(round(num_timesteps/num_frames))
    radius = 1.2*maximum([sqrt(p.x[1][1]^2 + p.x[1][2]^2) for p in planets])
    dot_scaling = 1/10
    lims = (-radius, radius)

    position_anim = @animate for t=1:gif_step:num_timesteps-1
        plot(title = "Solar System, t=$t", xlims=lims, ylims=lims)
        for p in planets
            x, y = p.x[t][1], p.x[t][2]
            scatter!([x], [y], markersize = log(p.m)*dot_scaling, label = p.name)
        end
    end
    gif(position_anim, "Week11/orbit.gif", fps = FPS)
end

function plot_velocities(planets)
    plot(title="velocities", xlabel="timestep", ylabel="abs. vel. [m/s]")
    for p in planets
        plot!(norm.(p.v), label=p.name)
    end
    png("Week11/velocities.png")
end

function plot_energies(p::Planet, sun::Planet)
    #=e_kin_planet(p) = [1/2 * p.m * norm.(p.v).^2]
    e_kin_system = sum(e_kin_planet(p) for p in planets)
    plot(e_kin_system, title="E_kin", xlabel="timestep", ylabel="energy [J]")
    png("Week11/energy.png")
    =#
    e_kin = 1/2 * p.m * norm.(p.v).^2
    e_pot = G * p.m * sun.m * norm.(p.x - sun.x).^(-1) # needs minus but idk...
    # TODO normalize
    e_tot = (e_kin + e_pot)

    plot(title="Energy $(p.name)", xlabel="timestep", ylabel="energy (normalized)")
    plot!(e_kin, label="E_kin")
    plot!(e_pot, label="E_pot")
    plot!(e_tot, label="e_tot")
    png("Week11/energy_$(p.name).png")
end

function main()
    # SIMULATION PARAMETERS
    dt = 5 * 60 # [s] timestep
    one_year = 12*30*24*60*60 # [s]
    T = 4 * one_year # [s] end time
    num_timesteps = Int(round(T/dt))
    animation_duration = 3

    @show num_timesteps
    @show dt

    println("initializing solar system...")
    planets = initial_config(num_timesteps)
    names = [p.name for p in planets]
    get_planet = Dict(names .=> planets)

    println("running simulation...")
    integrate_leapfrog!(planets, num_timesteps, dt)
    println("finished simulation.")
    println("animating orbit...")
    animate_orbit(planets, animation_duration) 
    println("plotting earth energies...")
    plot_energies(get_planet["Earth"], get_planet["Sun"]) #TODO hardcoded and fragile. not guud.
    println("plotting velocities...")
    plot_velocities(planets)
    println("done.")  
end

main() 
