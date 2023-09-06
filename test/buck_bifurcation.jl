include("../src/factorio.jl")

# E_range = 15:40
E_range = 15:0.001:40
schedule = DataFrame(idx = eachindex(E_range), E = E_range)

xdots = Float64[]; ydots = Float64[]
for dr in ProgressBar(eachrow(schedule))
    data = factory_buck(dr.idx, dr.E)
    
    idx_sampled = diff(data.Vr) .< 0
    sampledV = data[Not(1), :V][idx_sampled]
    append!(xdots, fill(dr.E, length(sampledV)))
    append!(ydots, sampledV)
end
@time a1 = scatter(xdots, ydots,
xlabel = L"E", ylabel = L"V",
label = :none, msw = 0, color = :black, ms = 0.5, alpha = 0.5, size = (700, 480));

png(a1, "buck_bifurcation")

# E_range = [20, 30.5, 31.5, 32.1, 33, 40]
# schedule = DataFrame(idx = eachindex(E_range), E = E_range)
# for dr in ProgressBar(eachrow(schedule))
#     data = factory_buck(dr.idx, dr.E)
# end


# schedule[100:800,:]
# dr = eachrow(schedule)[92]
# # @time data = factory_buck(dr.idx, dr.E); # to check compilation time
# @time data = factory_buck(dr.idx, dr.E)
# plot(data.V, data.I, alpha = 0.9, legend = :best, label = :none, xlabel = L"V", ylabel = L"I", color = :black, size = (400, 400), title = L"E = 30.5")
# png("a2")

# idx_sampled = diff(data.Vr) .< 0
# sampledV = data[Not(1), :V][idx_sampled]
# scatter!(data[Not(1), :V][idx_sampled], data[Not(1), :I][idx_sampled], alpha = 0.9, label = "Sampled", color = :black, ms = 1, msw = 0)
# # idx_jump = (diff(data.dI)) .< -100;
# # scatter!(data.V[Not(end)][idx_jump], data.I[Not(end)][idx_jump], alpha = 0.9)