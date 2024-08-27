include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

schedules = CSV.read("bifurcation/buck_schedules.csv", DataFrame)[1:1:end,:]
function buck_recovery()
    vrbl = [:dV, :dI], [:V, :I]
    cnfg = (; N = 1)
    dt = 1e-7
    θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;

    @threads for dr in eachrow(schedules) # dr = eachrow(schedules)[end]
        filename = "bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv"
        isfile(filename) && continue
        data = CSV.read(replace(filename, "buck_rcvd" => "buck"), DataFrame);

        add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
        Dtree = dryad(data, [:V, :I, :Vr])
        # Dtree = dryad(data, last(vrbl))

        ic = collect(data[1, last(vrbl)])
        try
            ŷ = DataFrame(solve(f_, ic, dt, data.Vr, Dtree, data.Vr), last(vrbl))
            CSV.write(filename, ŷ)
        catch
            @error "Error: $(lpad(dr.idx, 5, '0'))"
        end
    end
end
buck_recovery()

_Vr = CSV.read("bifurcation/buck/00001.csv", DataFrame).Vr
idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedules)
    filename = "bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) || continue
    data = CSV.read(filename, DataFrame)

    idx_sampled = diff(_Vr) .< 0
    sampledV = data[Not(1), :V][idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledV)))
    append!(hrzn, fill(dr.E, length(sampledV)))
    append!(vrtc, sampledV)
end
scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1)
CSV.write("bifurcation/buck_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
png("bifurcation/buck_recovered.png")

# schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)
# function soft_recovery()
#     vrbl = [:dt, :du, :dv], [:t, :u, :v]
#     cnfg = (; f_ = [cospi, sign], λ = 1e-2)
#     # cnfg = (; M = 3, λ = 1e-2)
#     dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

#     @threads for dr in eachrow(schedules)[1:1:end,:] # dr = eachrow(schedules)[192]
#         filename = "bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#         isfile(filename) && continue
#         data = CSV.read(replace(filename, "soft_rcvd" => "soft"), DataFrame);

#         add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 15~
#         f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
#         Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
        
#         try
#             ic = collect(data[1, last(vrbl)])
#             ŷ = DataFrame(solve(f_, ic, dt, data.t, Dtree), last(vrbl))
#             CSV.write(filename, ŷ)
#             CSV.write("G:/temp.csv", ŷ)
#         catch
#             @error "Error: $(lpad(dr.idx, 5, '0'))"
#         end
#     end
# end
# soft_recovery()

# idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
# @time for dr in eachrow(schedules)[1:1:end]
#     try
#         filename = "bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#         isfile(filename) || continue
#         data = CSV.read(filename, DataFrame)

#         idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
#         sampledv = data[Not(1), :v][idx_sampled]
#         append!(idcs, fill(dr.idx, length(sampledv)))
#         append!(hrzn, fill(dr.d, length(sampledv)))
#         append!(vrtc, sampledv)
#         print(dr.idx, "-")
#     catch
#         @error "Error: $(lpad(dr.idx, 5, '0'))"
#     end
# end
# CSV.write("bifurcation/soft_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = .1);
# png("bifurcation/soft_recovered.png")
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = .1, ylims = [-1, 1]);
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [430, 440])
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [440, 450], ticks = 440:450)
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [550, 570])
# temp = DataFrame(; idcs, vrtc, hrzn)
# temp = temp[temp.idcs .∉ Ref([439, 448, 560]), :]
# scatter(temp.hrzn, temp.vrtc, ms = 1, legend = :none, msw = 0, ma = 1);
# png("temp")
# @info ""

# schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)[1:1:end,:]
# function hrnm_recovery()
#     vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
#     cnfg = (; N = 3, f_ = [cos])
#     dt = 1e-3

#     @threads for dr in eachrow(schedules) # dr = eachrow(schedules)[171]
#         filename = "bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#         isfile(filename) && continue
#         data = CSV.read(replace(filename, "hrnm_rcvd" => "hrnm"), DataFrame);

#         add_subsystem!(data, vrbl, cnfg; θ1 = 3e-2, θ2 = 1e-10, min_rank = 20); # 30 sec
#         f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
#         Dtree = dryad(data, last(vrbl))

#         ic = collect(data[1, last(vrbl)])
#         try
#             ŷ = DataFrame(solve(f_, ic, dt, data.t, Dtree), last(vrbl))
#             CSV.write(filename, ŷ)
#         catch
#             @error "Error: $(lpad(dr.idx, 5, '0'))"
#         end
#     end
# end
# hrnm_recovery()

# idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
# @showprogress for dr in eachrow(schedules)
#     filename = "bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#     isfile(filename) || continue
#     data = CSV.read(filename, DataFrame)

#     idx_sampled = abs.(diff(diff(data.z) ./ 1e-3)) .> 0.1
#     sampledx = data[Not(1, end), :x][idx_sampled]
#     append!(idcs, fill(dr.idx, length(sampledx)))
#     append!(hrzn, fill(dr.bp, length(sampledx)))
#     append!(vrtc, sampledx)
# end
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1)
# CSV.write("bifurcation/hrnm_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# png("bifurcation/hrnm_recovered.png")