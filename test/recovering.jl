include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

schedule = CSV.read("G:/DDM/bifurcation/buck_schedule.csv", DataFrame)[1:end,:]

@showprogress for dr in eachrow(schedule)
    isfile("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv") && continue

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

    CSV.write("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv", ŷ)
    # idx_sampled = diff(data.Vr) .< 0
    # sampledV = ŷ[Not(1), :V][idx_sampled]
    # append!(vrtc, fill(dr.E, length(sampledV)))
    # append!(hrzn, sampledV)
    # CSV.write("G:/DDM/bifurcation/buck_recovered.csv", DataFrame(; vrtc, hrzn))
end

_Vr = CSV.read("G:/DDM/bifurcation/buck/00001.csv", DataFrame).Vr
vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedule)
    data = CSV.read("G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    idx_sampled = diff(_Vr) .< 0
    sampledV = data[Not(1), :V][idx_sampled]
    append!(vrtc, fill(dr.E, length(sampledV)))
    append!(hrzn, sampledV)
end
scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1);
png("G:/DDM/bifurcation/buck_recovered.png")


schedule = CSV.read("G:/DDM/bifurcation/soft_schedule.csv", DataFrame)[1:1:end,:]

@showprogress for dr in eachrow(schedule)
    isfile("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv") && continue

    data = CSV.read("G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    f = SINDy(data[rand(1:nrow(data), 22), :], [:dt, :du, :dv], [:t, :u, :v], verbose=false)
    while true
        sampled = rand(1:nrow(data), 22)
        f = SINDy(data[sampled, :], [:dt, :du, :dv], [:t, :u, :v], M = 3, verbose=false)
        if f.MSE < 1e-10
            print(f, ["t", "u", "v"])
            break
        end
    end

    error_ = norm.(eachrow(Matrix(data[:, [:dt, :du, :dv]])) .- f.(eachrow(Matrix(data[:, [:t, :u, :v]]))))
    bit_alien = error_ .> 1e-4

    subsystem = ones(Int64, nrow(data));
    subsystem[bit_alien] .= 2;
    data[:, :subsystem] = subsystem;

    test = data
    _trng = data
    gdf_ = groupby(_trng, :subsystem)
    f_ = [SINDy(gdf, [:dt, :du, :dv], [:t, :u, :v], M = 3) for gdf in gdf_]
    
    labels = _trng.subsystem
    features = Matrix(_trng[:, [:t, :u, :v]])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
    println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

    dt = 1e-5
    x = collect(test[1, [:t, :u, :v]])
    y = test

    try
        ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "u", "v"]);
        CSV.write("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv", ŷ)
    catch
        println("Error: $(lpad(dr.idx, 5, '0'))")
    end
end

idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedule)
    data = CSV.read("G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
    sampledv = data[Not(1), :v][idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledv)))
    append!(vrtc, fill(dr.d, length(sampledv)))
    append!(hrzn, sampledv)
end
CSV.write("G:/DDM/bifurcation/soft_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
scatter(vrtc, hrzn, ms = 1, legend = :none, msw = 0, ma = 0.1, ylims = [-1, 1]);
png("G:/DDM/bifurcation/soft_recovered.png")
# rm.(string.("G:/DDM/bifurcation/soft_rcvd/", lpad.(unique(idcs[abs.(hrzn) .> .64]), 5, '0'), ".csv"))

schedule = CSV.read("G:/DDM/bifurcation/hrnm_schedule.csv", DataFrame)[1:1:end,:]

dr = first(eachrow(schedule))
theset = [1,5]
plot(data.z)
plot(data.dz)
plot(normeddf)
# @showprogress for dr in eachrow(schedule)
    isfile("G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv") && continue

    data = CSV.read("G:/DDM/bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv", DataFrame)

    normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, [:dx, :dy, :dz]]), dims = 1)))
    jumpt = [1; findall(normeddf .> 2e-2)]
    # plot(trng.t, trng.z)
    # hline!([1, -1], color = :blue)
    # vline!(trng.t[jumpt], color = :red, ls = :dash)
    
    subsystem = zeros(Int64, nrow(data));
    sets = collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1))); pop!(sets); sets = UnitRange.(first.(sets), last.(sets))
    for id_subsys in 1:4
        candy_ = []
        for n in 2:3
            for theset = combinations(eachindex(sets), n)
                # if idx_long ∉ theset continue end
                if 1 ∈ diff(theset) continue end
                print(theset)
    
                sliced = data[reduce(vcat, sets[theset]), :]
                X = Θ(sliced[:, [:t, :x, :y, :z]], N = 3, f_ = [cos])
                if rank(X) ≥ size(X, 2)
                    candy = SINDy(sliced,
                            [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
                            N = 3, f_ = [cos])
                    print(candy, ["t", "x", "y", "z"])
                    push!(candy_, theset => candy.MSE)
                    if candy.MSE < 1e-28 break end
                end
            end
        end
        picked = first(candy_[argmin(last.(candy_))])
        f = SINDy(data[reduce(vcat, sets[picked]), :],
            [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
            N = 3, f_ = [cos])
        print(f, ["t", "x", "y", "z"])
    
        @time error_ = norm.(eachrow(Matrix(data[:, [:dt, :dx, :dy, :dz]])) .- f.(eachrow(Matrix(data[:, [:t, :x, :y, :z]]))))
        bit_alien = error_ .> 1e-8
        # scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")
        idx_alien = findall(.!bit_alien)
        subsystem[idx_alien] .= id_subsys
        # plot(data.u[1:100:100000], data.v[1:100:100000], color=data.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))
    
        # sets = filter(!isempty, [setdiff(set, idx_alien) for set in sets])
        # sets = filter(x -> isdisjoint(x, idx_alien), sets)
        idx_unlabled = []
        for i in eachindex(sets)
            if first(sets[i]) ∉ idx_alien
                push!(idx_unlabled, i)
            end
        end
        sets = sets[idx_unlabled]
        if isempty(sets)
            @info "all data points are exhausted!"
            break
        end
    end
    data[:, :subsystem] = subsystem;
    

    test = data
    _trng = data
    gdf_ = groupby(_trng, :subsystem)
    f_ = [SINDy(gdf, [:dt, :dx, :dy, :dz], [:t, :x, :y, :z], N = 3, f_ = [cos]) for gdf in gdf_]

    labels = _trng.subsystem
    features = Matrix(_trng[:, [:t, :x, :y, :z]])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
    println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

    dt = 1e-3
    x = collect(test[1, [:t, :x, :y, :z]])
    y = test
    ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "x", "y", "z"])

    CSV.write("G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv", ŷ)
# end

vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedule)
    idx_sampled = abs.(diff(data.dz)) .> 0.1
    sampledx = data[Not(1), :x][idx_sampled]
    append!(vrtc, fill(dr.f, length(sampledx)))
    append!(hrzn, sampledx)
end
scatter(vrtc, hrzn, color = :black, ms = 1, legend = :none, msw = 0, ma = 1)
png("G:/DDM/bifurcation/hrnm_recovered.png")