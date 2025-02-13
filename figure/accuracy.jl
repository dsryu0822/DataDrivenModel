

rfrc = CSV.read("result/hrnm_lyapunov.csv", DataFrame)
rcvd = CSV.read("result/hrnm_lyapunov_rcvd.csv", DataFrame)
rfrc = Matrix(rfrc[:, 3:5])
rcvd = Matrix(rcvd[:, 3:5])

tol_ = logrange(1e-6, 7e-2, 100)
acc_ = []
for tol = tol_
    _rfrc = sign.(rfrc) .* (abs.(rfrc) .> tol)
    _rcvd = sign.(rcvd) .* (abs.(rcvd) .> tol)
    push!(acc_, sum(_rfrc .== _rcvd) / lastindex(_rfrc))
end
acc_
plot(acc_)
maximum(acc_)

[0.9355322338830585, 0.9943820224719101, 0.997002997002997]