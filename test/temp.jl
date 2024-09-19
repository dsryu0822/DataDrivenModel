include("../core/header.jl")
using LaTeXStrings

data = CSV.read("lyapunov/soft_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/soft_lyapunov_rcvd.csv", DataFrame)
plot(xlabel = L"d", ylabel = L"\lambda", legend = :none, title = "Lyapunov spectrum of soft impact model")
plot!(data.bp, data.λ1, color = :gray)
plot!(data.bp, data.λ2, color = :gray)
plot!(data.bp, data.λ3, color = :gray)
plot!(rcvd.bp, rcvd.λ1, color = :red)
plot!(rcvd.bp, rcvd.λ2, color = :red)
plot!(rcvd.bp, rcvd.λ3, color = :red)
png("lyapunov/soft_lyapunov")


data = CSV.read("lyapunov/hrnm_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/hrnm_lyapunov_rcvd.csv", DataFrame)
plot(xlabel = L"f", ylabel = L"\lambda", legend = :none, title = "Lyapunov spectrum of HR model")
plot!(data.bp, data.λ1, color = :gray)
plot!(data.bp, data.λ2, color = :gray)
plot!(data.bp, data.λ3, color = :gray)
plot!(rcvd.bp, rcvd.λ1, color = :red)
plot!(rcvd.bp, rcvd.λ2, color = :red)
plot!(rcvd.bp, rcvd.λ3, color = :red)
png("lyapunov/hrnm_lyapunov")


data = CSV.read("lyapunov/gear_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/gear_lyapunov_rcvd.csv", DataFrame)
plot(xlabel = L"Fe", ylabel = L"\lambda", legend = :none, title = "Lyapunov spectrum of gear model")
plot!(data.bp, data.λ1, color = :gray)
plot!(data.bp, data.λ2, color = :gray)
plot!(data.bp, data.λ3, color = :gray)
plot!(rcvd.bp, rcvd.λ1, color = :red)
plot!(rcvd.bp, rcvd.λ2, color = :red)
plot!(rcvd.bp, rcvd.λ3, color = :red)
png("lyapunov/gear_lyapunov")

bfcn = CSV.read("bifurcation/gear_bifurcation_rcvd.csv", DataFrame)
scatter(bfcn.hrzn, bfcn.vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, color = :black);
png("bifurcation/gear_bifurcation_rcvd.png")