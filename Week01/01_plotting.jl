
using Plots
using LinearAlgebra
using Images
using LsqFit

x = LinRange(0, 2π, 100)

# The dot indicates that the function is applied elementwise, referred to as "broadcasting".
# The semicolon ';' at the end prevents the notebook to print the return value.
y = sin.(x)
plot(x, y, label = "sin(x)")
y = exp.(x) / 300
plot!(x, y, label = "C*eˣ")
y = cos.(x)
plot!(x, y, label = "cos(x)")

plot!(title = "A plot title", xlabel = "x", ylabel = "y")
png("plot1")

