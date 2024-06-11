include("../core/DDM.jl")
include("../core/factorio.jl")
using ProgressMeter
packages = [:CSV, :DataFrames, :Plots, :Colors, :ColorSchemes]
for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm
using Base.Threads: @threads


# E_range = 15:0.01:40
# schedules = DataFrame(idx = eachindex(E_range), E = E_range)
# vrtc = Float64[]; hrzn = Float64[]
# CSV.write("G:/DDM/bifurcation/buck_schedules.csv", schedules)
# @showprogress @threads for dr in eachrow(schedules)
#     data = factory_buck(DataFrame, dr.idx, dr.E, tspan = [0, 0.30])
#     data = data[29(nrow(data) ÷ 30):end, :]
#     CSV.write("G:/DDM/bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv", data)
    
#     idx_sampled = diff(data.Vr) .< 0
#     sampledV = data[Not(1), :V][idx_sampled]
#     append!(hrzn, fill(dr.E, length(sampledV)))
#     append!(vrtc, sampledV)
# end
# CSV.write("G:/DDM/bifurcation/buck_bifurcation.csv", DataFrame(; vrtc, hrzn))
# scatter(hrzn, vrtc, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
# png("G:/DDM/bifurcation/buck_bifurcation.png")

d_range = 0.1:0.0001:0.3
schedules = DataFrame(idx = eachindex(d_range), d = d_range)
vrtc = Float64[]; hrzn = Float64[]
CSV.write("G:/DDM/lyapunov/soft_schedules.csv", schedules)
@showprogress @threads for dr in eachrow(schedules)
    # 30~50: enough to get bifurcation
    # 30~100: not enough to get lyaunov exponent
    data = factory_soft(DataFrame, dr.idx, dr.d, tspan = [0, 150])
    data = data[30(nrow(data) ÷ 150):end, :]
    CSV.write("G:/DDM/lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv", data)
    
    # idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
    # sampledv = data[Not(1), :v][idx_sampled]
    # append!(hrzn, fill(dr.d, length(sampledv)))
    # append!(vrtc, sampledv)
end
# CSV.write("G:/DDM/bifurcation/soft_bifurcation.csv", DataFrame(; vrtc, hrzn))
# scatter(hrzn, vrtc, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
# png("G:/DDM/bifurcation/soft_bifurcation.png")

# f_range = 0:0.001:1
# schedules = DataFrame(idx = eachindex(f_range), f = f_range)
# vrtc = Float64[]; hrzn = Float64[]
# CSV.write("G:/DDM/bifurcation/hrnm_schedules.csv", schedules)
# @showprogress @threads for dr in eachrow(schedules)
#     data = factory_hrnm(DataFrame, dr.idx, dr.f, tspan = [0, 1500])
#     data = data[1000(nrow(data) ÷ 1500):end, :]
#     CSV.write("G:/DDM/bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv", data)

#     idx_sampled = abs.(diff(data.dz)) .> 0.1
#     sampledx = data[Not(1), :x][idx_sampled]
#     append!(hrzn, fill(dr.f, length(sampledx)))
#     append!(vrtc, sampledx)
# end
# CSV.write("G:/DDM/bifurcation/hrnm_bifurcation.csv", DataFrame(; vrtc, hrzn))
# scatter(hrzn, vrtc, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
# png("G:/DDM/bifurcation/hrnm_bifurcation.png")