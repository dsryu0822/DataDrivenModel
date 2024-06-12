using DecisionTree, Random, StatsBase, Dates; @info now()
using Base.Threads: @threads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")

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
# vrtc = Float64[]; hrzn = Float64[]
# CSV.write("G:/DDM/bifurcation/soft_schedules.csv", schedules)
# @showprogress @threads for dr in eachrow(schedules)
#     # 30~50: enough to get bifurcation
#     # 30~100: not enough to get lyaunov exponent
#     data = factory_soft(DataFrame, dr.idx, dr.d, tspan = [0, 50])
#     data = data[30(nrow(data) ÷ 50):end, :]
#     CSV.write("G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv", data)
    
#     # idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
#     # sampledv = data[Not(1), :v][idx_sampled]
#     # append!(hrzn, fill(dr.d, length(sampledv)))
#     # append!(vrtc, sampledv)
# end
# CSV.write("G:/DDM/bifurcation/soft_bifurcation.csv", DataFrame(; vrtc, hrzn))
# scatter(hrzn, vrtc, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
# png("G:/DDM/bifurcation/soft_bifurcation.png")

vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2)
dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
# lastpoint = []
@showprogress @threads for dr = eachrow(schedules)
    filename = "G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
    !isfile(filename) && continue

    data = CSV.read(filename, DataFrame)
    # push!(lastpoint, collect(data[end, 1:3]))
    
    if !hasproperty(data, :subsystem)
        add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
        CSV.write(filename, data)
    end
end
# CSV.write("G:/DDM/lyapunov/soft_schedules_cache.csv", [schedules DataFrame(stack(lastpoint)', [:t, :u, :v])])

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