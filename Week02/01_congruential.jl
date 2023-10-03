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
    y = copy(x)
    popfirst!(y)
    pop!(x)
    scatter(x, y, label = "congruential rn")
    plot!(title = "Square Test", xlabel = "x_i", ylabel = "x_(i+1)")
    png("Week02/square_test")
end

function cube_test(x)
    y = copy(x)
    popfirst!(y)
    z = copy(y)
    popfirst!(z)
    pop!(x); pop!(x)
    pop!(y)
    if size(x) == size(y) == size(z)
        println("all vectors same size")        
    end
    scatter(x, y, z, label = "congruential rn")
    plot!(title = "Cube Test", xlabel = "x_i", ylabel = "x_(i+1)", zlabel = "x_(i+2)")
    png("Week02/cube_test")
end 

n = 1000
square_test(congruential_rn(n)) # create one extra for i+1 comparison
cube_test(congruential_rn(n))
