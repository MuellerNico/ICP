using LsqFit
using Plots
using LinearAlgebra

xdata = [ 15.2; 19.9; 2.2; 11.8; 12.1; 18.1; 11.8; 13.4; 11.5; 0.5;
18.0; 10.2; 10.6; 13.8; 4.6; 3.8; 15.1; 15.1; 11.7; 4.2 ]
ydata = [ 0.73; 0.19; 1.54; 2.08; 0.84; 0.42; 1.77; 0.86; 1.95; 0.27;
0.39; 1.39; 1.25; 0.76; 1.99; 1.53; 0.86; 0.52; 1.54; 1.05 ]

beta0 = [1.0, 1.0, 1.0]

@. f(x, beta) = beta[1] * (x./beta[2]) * exp(- (x/beta[2]) ^ (beta[3]) )

fit = curve_fit(f, xdata, ydata, beta0, lower=[-Inf, 0, -Inf])
beta = fit.param
err = stderror(fit)

#beta1 = coef(fit)
println("beta = $beta")

xs = LinRange(0, 20, 100)
ys = f(xs, beta)

scatter(xdata, ydata, label = "raw data")
plot!(xs, ys, label = "fitting function")
png("fitting")

