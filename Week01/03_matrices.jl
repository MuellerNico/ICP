using LinearAlgebra

function randmatrix(n)
    #A = [rand(Float32) for i=1:n, j=1:n]
    A = rand(n,n)
    Q,R = qr(A)
    return Q
end

function matgen(n, c)
    U = randmatrix(n)
    V = randmatrix(n)
    s = LinRange(1/c, 1, n)
    σ = Diagonal(s)
    return U*σ*V
end 

function is_orthogonal(M)
    return maximum(abs.(M*M' - I)) < 1e-14
end

c = 15
n = 4
A = matgen(n, c)
println("Matrix is orthogonal: $(is_orthogonal(A))")
show(stdout, "text/plain", A)
println()
k = norm(A) * norm(inv(A))
println("k calc = $k, k spec = $c")