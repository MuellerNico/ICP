using Distributed
using Dates
using DelimitedFiles
using JSON
@everywhere include("shallow_water.jl")


function benchmark_parallel(L)
    @everywhere g = 9.81
    @everywhere b = 0.3
    @everywhere ν = 0.3

    @everywhere Δx = 0.5
    @everywhere Δt = 0.01

    @everywhere A = 1.
    @everywhere σ = 5.

    @everywhere n_steps = 200;

    @everywhere process_topologies = Dict(
        1 => (1, 1),
        2 => (1, 2),
        3 => (1, 3),
        4 => (1, 4),
        5 => (1, 5),
        6 => (1, 6),
    )

    @everywhere n_proc_rows, n_proc_cols = process_topologies[nprocs()]
    @everywhere neighbours = get_neighbours(n_proc_rows, n_proc_cols);

    # Copy the function argument to every process.
    # Taken from: https://stackoverflow.com/a/37027837
    @eval @everywhere L=$L

    # Total linear grid size; each process only has a part of this!
    @everywhere N = Int64(ceil(L / Δx))

    @show L

    t_start = Dates.value(Dates.now())
    @everywhere grid = simulate_parallel(
        N, L, Δx, Δt, g, b, ν, n_steps, A, σ, n_proc_rows, n_proc_cols
    );
    t_end = Dates.value(Dates.now())

    return t_end - t_start
end

L_values = [60, 120, 240, 480]

timings = Dict()
for L in L_values
    timings[L] = benchmark_parallel(L)
end

println(timings)
open("benchmark_parallel_nproc_$(nprocs()).json", "w") do file
    write(file, JSON.json(timings))
end
