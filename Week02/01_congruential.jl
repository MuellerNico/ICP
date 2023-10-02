using Plots

function congruential_rn(n)
    c = 3
    p = 31
    x = 4 #x1
    rn = zeros(n)
    for i=1:n
        x = c*x % p
        rn[i] = x
    end
    return rn
end

function square_test(x)
    scatter(x, popfirst!(x), label = "congruential rn") # TODO ....
    plot!(title = "Square Test", xlabel = "x_i", ylabel = "x_(i+1)")
    png("Week02/square_test")
end

function cube_test(x)
end 

n = 100
square_test(congruential_rn(n)) # create one extra for i+1 comparison
cube_test(congruential_rn(n))
