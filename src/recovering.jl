include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random
using Base.Threads: @threads # Base.Threads.nthreads()

# plot(data[sets[4], :z])
# candy = SINDy(data[meaningful,:], vrbl...; cnfg...); print(candy, last(vrbl))

candy_ = [SINDy(data[a_set,:], vrbl...; cnfg...) for a_set in sets]
print.(candy_, Ref(last(vrbl)))
[rank(Θ(Matrix(data[a_set,last(vrbl)]); cnfg...)) for a_set in sets]

Θ(Matrix(data[subsystem .== 0,last(vrbl)]); cnfg...)
rank(Θ(Matrix(data[subsystem .== 0,last(vrbl)]); cnfg...))
data[subsystem .== 0, :]
plot(data[meaningful,:z])
plot(data[sets[argmax(rank_)],:z])

function add_subsystem!(data, vrbl, cnfg; θ1 = 2e-2, θ2 = 1e-10, min_rank = 0)
    normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1))) # plot(normeddf)
    jumpt = [1; findall(normeddf .> θ1)]
    sets = filter(!isempty, collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1)))); pop!(sets); sets = UnitRange.(first.(sets), last.(sets))

    subsystem = zeros(Int64, nrow(data));
    for id_subsys = 1:4; # id_subsys = 1; id_subsys = 2; id_subsys = 3
        if sets |> isempty
            # @info "all data points are labeled!"
            break
        end
        rank_ = [rank(Θ(Matrix(data[a_set,last(vrbl)]); cnfg...)) for a_set in sets]
        if maximum(rank_) < min_rank
            candy = SINDy(data[iszero.(subsystem),:], vrbl...; cnfg...); #print(candy, last(vrbl))
        else
            meaningful = sets[argmax(rank_)]
            candy = SINDy(data[meaningful,:], vrbl...; cnfg...); #print(candy, last(vrbl))
        end

        idx_blank = findall(iszero.(subsystem))
        idx_blank = idx_blank[norm.(eachrow(Matrix(data[idx_blank, first(vrbl)])) .- candy.(eachrow(Matrix(data[idx_blank, last(vrbl)])))) .< θ2]
        subsystem[idx_blank] .= id_subsys
        sets = sets[rand.(sets) .∉ Ref(idx_blank)]
    end
    data[!, :subsystem] = subsystem;
    return data
end
function dryad(data, vrbl) # fairy of tree and forest
    labels = data.subsystem
    features = Matrix(data[:, last(vrbl)])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(data.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(data.subsystem, features, rng = argmax(acc_))
    # println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")
    return Dtree
end

# schedule = CSV.read("G:/DDM/bifurcation/buck_schedule.csv", DataFrame)[1:end,:]

# @showprogress for dr in eachrow(schedule)
#     isfile("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv") && continue

#     data = CSV.read("G:/DDM/bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
#     f = SINDy(data[rand(1:nrow(data), 4), :], [:dV, :dI], [:V, :I], verbose=false)
#     while true
#         sampled = rand(1:nrow(data), 4)
#         f = SINDy(data[sampled, :], [:dV, :dI], [:V, :I], verbose=false)
#         if f.MSE < 1e-10
#             # print(f, ["V", "I"])
#             break
#         end
#     end

#     error_ = norm.(eachrow(Matrix(data[:, [:dV, :dI]])) .- f.(eachrow(Matrix(data[:, [:V, :I]]))))
#     bit_alien = error_ .> 1e-4

#     subsystem = ones(Int64, nrow(data));
#     subsystem[bit_alien] .= 2;
#     data[:, :subsystem] = subsystem;

#     test = data
#     _trng = data
#     gdf_ = groupby(_trng, :subsystem)
#     f_ = [SINDy(gdf, [:dV, :dI], [:V, :I]) for gdf in gdf_]

#     labels = _trng.subsystem
#     features = Matrix(_trng[:, [:V, :I, :Vr]])
#     acc_ = []
#     for seed in 1:10
#         Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
#         push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
#         if maximum(acc_) ≈ 1 break; else print("█") end
#     end
#     Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
#     # println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

#     dt = 1e-7
#     x = collect(test[1, [:V, :I]])
#     y = test
#     ŷ = DataFrame(solve(f_, x, dt, y.Vr, Dtree, y.Vr), ["V", "I"]);

#     CSV.write("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv", ŷ)
#     # idx_sampled = diff(data.Vr) .< 0
#     # sampledV = ŷ[Not(1), :V][idx_sampled]
#     # append!(vrtc, fill(dr.E, length(sampledV)))
#     # append!(hrzn, sampledV)
#     # CSV.write("G:/DDM/bifurcation/buck_recovered.csv", DataFrame(; vrtc, hrzn))
# end

# _Vr = CSV.read("G:/DDM/bifurcation/buck/00001.csv", DataFrame).Vr
# vrtc = Float64[]; hrzn = Float64[]
# @showprogress for dr in eachrow(schedule)
#     data = CSV.read("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
#     idx_sampled = diff(_Vr) .< 0
#     sampledV = data[Not(1), :V][idx_sampled]
#     append!(vrtc, fill(dr.E, length(sampledV)))
#     append!(hrzn, sampledV)
# end
# scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
# png("G:/DDM/bifurcation/buck_recovered.png")


# schedule = CSV.read("G:/DDM/bifurcation/soft_schedule.csv", DataFrame)[1:1:end,:]

# @showprogress for dr in eachrow(schedule)
#     isfile("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv") && continue

#     data = CSV.read("G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
#     f = SINDy(data[rand(1:nrow(data), 22), :], [:dt, :du, :dv], [:t, :u, :v], verbose=false)
#     while true
#         sampled = rand(1:nrow(data), 22)
#         f = SINDy(data[sampled, :], [:dt, :du, :dv], [:t, :u, :v], M = 3, verbose=false)
#         if f.MSE < 1e-10
#             print(f, ["t", "u", "v"])
#             break
#         end
#     end

#     error_ = norm.(eachrow(Matrix(data[:, [:dt, :du, :dv]])) .- f.(eachrow(Matrix(data[:, [:t, :u, :v]]))))
#     bit_alien = error_ .> 1e-4

#     subsystem = ones(Int64, nrow(data));
#     subsystem[bit_alien] .= 2;
#     data[:, :subsystem] = subsystem;

#     test = data
#     _trng = data
#     gdf_ = groupby(_trng, :subsystem)
#     f_ = [SINDy(gdf, [:dt, :du, :dv], [:t, :u, :v], M = 3) for gdf in gdf_]
    
#     labels = _trng.subsystem
#     features = Matrix(_trng[:, [:t, :u, :v]])
#     acc_ = []
#     for seed in 1:10
#         Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
#         push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
#         if maximum(acc_) ≈ 1 break; else print("█") end
#     end
#     Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
#     println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

#     dt = 1e-5
#     x = collect(test[1, [:t, :u, :v]])
#     y = test

#     try
#         ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "u", "v"]);
#         CSV.write("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv", ŷ)
#     catch
#         println("Error: $(lpad(dr.idx, 5, '0'))")
#     end
# end

# idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
# @showprogress for dr in eachrow(schedule)
#     data = CSV.read("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
#     idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
#     sampledv = data[Not(1), :v][idx_sampled]
#     append!(idcs, fill(dr.idx, length(sampledv)))
#     append!(hrzn, fill(dr.d, length(sampledv)))
#     append!(vrtc, sampledv)
# end
# CSV.write("G:/DDM/bifurcation/soft_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, ylims = [-1, 1]);
# png("G:/DDM/bifurcation/soft_recovered.png")
# # rm.(string.("G:/DDM/bifurcation/soft_rcvd/", lpad.(unique(idcs[abs.(hrzn) .> .64]), 5, '0'), ".csv"))

schedule = CSV.read("G:/DDM/bifurcation/hrnm_schedule.csv", DataFrame)[1:1:end,:]
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3

@showprogress for dr in eachrow(schedule) # dr = eachrow(schedule)[1]
    filename = "G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) && continue
    data = CSV.read(replace(filename, "hrnm_rcvd" => "hrnm"), DataFrame);

    @time add_subsystem!(data, vrbl, cnfg; θ1 = 3e-2, θ2 = 1e-5, min_rank = 20); # 30 sec
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
    Dtree = dryad(data, vrbl)

    ic = collect(data[1, last(vrbl)])
    try
        ŷ = DataFrame(solve(f_, ic, dt, data.t, Dtree), last(vrbl))
        CSV.write(filename, ŷ)
    catch
        @error "Error: $(lpad(dr.idx, 5, '0'))"
    end
end

idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedule)
    dr.idx == 56 && continue
    filename = "G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) || continue
    data = CSV.read(filename, DataFrame)

    idx_sampled = abs.(diff(diff(data.z) ./ 1e-3)) .> 0.1
    sampledx = data[Not(1, end), :x][idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledx)))
    append!(hrzn, fill(dr.f, length(sampledx)))
    append!(vrtc, sampledx)
end
CSV.write("G:/DDM/bifurcation/hrnm_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0)
scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, xlims = [1, 30], xticks = 1:60)
png("G:/DDM/bifurcation/hrnm_recovered.png")
