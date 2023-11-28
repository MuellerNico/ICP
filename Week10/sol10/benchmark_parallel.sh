# Serial version.
julia benchmark_parallel.jl

# Parallel versions.
for i in {1..5}
do
    julia -p "$i" benchmark_parallel.jl
done
python plot_benchmark.py
