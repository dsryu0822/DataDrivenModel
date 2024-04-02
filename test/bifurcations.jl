include("../src/factorio.jl")
using ProgressMeter
packages = [:CSV, :DataFrames, :Plots, :Colors, :ColorSchemes]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm

# E_range = 15:0.001:40
E_range = [collect(15:0.1:29); collect(29.001:0.001:40)]
schedule = DataFrame(idx = eachindex(E_range), E = E_range)

vrtc = Float64[]; hrzn = Float64[]
# CSV.write("G:/DDM/bifurcation/buck_schedule.csv", schedule)
@showprogress for dr in eachrow(schedule)
    data = factory_buck(DataFrame, dr.idx, dr.E, tspan = [0, 0.30])
    data = data[29(nrow(data) ÷ 30):end, :]
    # CSV.write("G:/DDM/bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv", data)
    
    idx_sampled = diff(data.Vr) .< 0
    sampledV = data[Not(1), :V][idx_sampled]
    append!(vrtc, fill(dr.E, length(sampledV)))
    append!(hrzn, sampledV)
end
# CSV.write("G:/DDM/bifurcation/buck_bifurcation.csv", DataFrame(; vrtc, hrzn))
# scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1)
# png("G:/DDM/bifurcation/buck_bifurcation.png")

# 
d_range = 0.1:0.001:0.3
schedule = DataFrame(idx = eachindex(d_range), d = d_range)

vrtc = Float64[]; hrzn = Float64[]
CSV.write("G:/DDM/bifurcation/soft_schedule.csv", schedule)
@showprogress for dr in eachrow(schedule)
    data = factory_soft(DataFrame, dr.idx, dr.d, tspan = [0, 50])
    data = data[45(nrow(data) ÷ 50):end, :]
    CSV.write("G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv", data)
    
    idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
    sampledu = data[Not(1), :v][idx_sampled]
    append!(vrtc, fill(dr.d, length(sampledu)))
    append!(hrzn, sampledu)
end
CSV.write("G:/DDM/bifurcation/soft_bifurcation.csv", DataFrame(; vrtc, hrzn))
scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1)
png("G:/DDM/bifurcation/soft_bifurcation.png")
