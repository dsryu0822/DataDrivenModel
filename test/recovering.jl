include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

schedule = CSV.read("G:/DDM/bifurcation/buck_schedule.csv", DataFrame)[1:10:end,:]

vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedule)
    data = CSV.read("G:/DDM/bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    f = SINDy(data[rand(1:nrow(data), 4), :], [:dV, :dI], [:V, :I], verbose=false)
    while true
        sampled = rand(1:nrow(data), 4)
        f = SINDy(data[sampled, :], [:dV, :dI], [:V, :I], verbose=false)
        if f.MSE < 1e-10
            # print(f, ["V", "I"])
            break
        end
    end

    error_ = norm.(eachrow(Matrix(data[:, [:dV, :dI]])) .- f.(eachrow(Matrix(data[:, [:V, :I]]))))
    bit_alien = error_ .> 1e-4

    subsystem = ones(Int64, nrow(data));
    subsystem[bit_alien] .= 2;
    data[:, :subsystem] = subsystem;

    test = data
    _trng = data
    gdf_ = groupby(_trng, :subsystem)
    f_ = [SINDy(gdf, [:dV, :dI], [:V, :I]) for gdf in gdf_]

    labels = _trng.subsystem
    features = Matrix(_trng[:, [:V, :I, :Vr]])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
    # println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

    dt = 1e-7
    x = collect(test[1, [:V, :I]])
    y = test
    ŷ = DataFrame(solve(f_, x, dt, y.Vr, Dtree, y.Vr), ["V", "I"]);

    idx_sampled = diff(data.Vr) .< 0
    sampledV = ŷ[Not(1), :V][idx_sampled]
    append!(vrtc, fill(dr.E, length(sampledV)))
    append!(hrzn, sampledV)
end
CSV.write("G:/DDM/bifurcation/buck_recovered.csv", DataFrame(; vrtc, hrzn))
scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1)
png("G:/DDM/bifurcation/buck_recovered.png")