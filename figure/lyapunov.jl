include("../core/header.jl")
using LaTeXStrings

## 원본은 ground thruth

data = CSV.read("lyapunov/soft_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/soft_lyapunov_rcvd.csv", DataFrame)
p3 = plot(legend = :none)
plot!(p3, data.bp, data.λ1, color = :gray)
plot!(p3, data.bp, data.λ2, color = :gray)
plot!(p3, data.bp, data.λ3, color = :gray)
plot!(p3, rcvd.bp, rcvd.λ1, color = 1)
plot!(p3, rcvd.bp, rcvd.λ2, color = 1)
plot!(p3, rcvd.bp, rcvd.λ3, color = 1)
# png("soft_lyapunov.png")


data = CSV.read("lyapunov/hrnm_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/hrnm_lyapunov_rcvd.csv", DataFrame)
plot(xlabel = L"f", ylabel = L"\lambda", legend = :none, title = "Lyapunov spectrum of HR model")
plot!(data.bp, data.λ1, color = :gray)
plot!(data.bp, data.λ2, color = :gray)
plot!(data.bp, data.λ3, color = :gray)
plot!(rcvd.bp, rcvd.λ1, color = 1)
plot!(rcvd.bp, rcvd.λ2, color = 1)
plot!(rcvd.bp, rcvd.λ3, color = 1)
png("lyapunov/hrnm_lyapunov")


data = CSV.read("lyapunov/gear_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/gear_lyapunov_rcvd.csv", DataFrame)
plot(xlabel = L"Fe", ylabel = L"\lambda", legend = :none, title = "Lyapunov spectrum of gear model")
plot!(data.bp, data.λ1, color = :gray)
plot!(data.bp, data.λ2, color = :gray)
plot!(data.bp, data.λ3, color = :gray)
plot!(rcvd.bp, rcvd.λ1, color = 1)
plot!(rcvd.bp, rcvd.λ2, color = 1)
plot!(rcvd.bp, rcvd.λ3, color = 1)
png("lyapunov/gear_lyapunov")
