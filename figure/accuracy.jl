

rfrc = CSV.read("result/hrnm_lyapunov.csv", DataFrame)
rcvd = CSV.read("result/hrnm_lyapunov_rcvd.csv", DataFrame)
rfrc = Matrix(rfrc[:, 3:5])
rcvd = Matrix(rcvd[:, 3:5])

# tol_ = logrange(1e-2, 1e-0, 100)
# tol_ = logrange(1e-2, 5e-2, 100)
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

[0.9881725803764785, 0.9883478984602581, 0.996003996003996]