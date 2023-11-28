M = zeros(3,3,3)
X = M
X[1,1,1] += 1
@show M

M = zeros(3,3,3)
X = @view M[1,:,:]
X[1,1] += 1
@show M
