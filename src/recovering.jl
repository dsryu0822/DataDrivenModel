include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, Clustering, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

function add_subsystem!(data, vrbl, cnfg; θ1 = 1e-1, θ2 = 1e-24, θ3 = 1e-10, min_rank = 0)
    normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    jumpt = [1; findall(normeddf .> θ1)]
    _sets = filter(!isempty, collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1)))); pop!(_sets); _sets = UnitRange.(first.(_sets), last.(_sets))
    sets = sort(union(_sets, UnitRange.(last.(_sets)[1:(end-1)] .+ 2, first.(_sets)[2:end] .- 2)))

    subsystem = zeros(Int64, nrow(data));
    for id_subsys = 1:4 # id_subsys = 1; id_subsys = 2; id_subsys = 3
        rank_ = [rank(Θ(Matrix(data[a_set,last(vrbl)]); cnfg...)) for a_set in sets]
        if maximum(rank_) < min_rank
            # candy = SINDy(data[iszero.(subsystem),:], vrbl...; cnfg...); print(candy, last(vrbl))
            for (A, B) = combinations(sets, 2)
                candy = SINDy([data[A, :]; data[B, :]], vrbl...; cnfg...)
                # A = first(sets); B = last(sets);
                if candy.MSE < θ2 print(candy, last(vrbl)); break else @info "move on to the next pair!" end
            end
        else
            meaningful = sets[argmax(rank_)]
            candy = SINDy(data[meaningful,:], vrbl...; cnfg...); print(candy, last(vrbl))
        end

        idx_blank = findall(iszero.(subsystem))
        residual = norm.(eachrow(Matrix(data[idx_blank, first(vrbl)])) .- candy.(eachrow(Matrix(data[idx_blank, last(vrbl)]))))
        # scatter(residual[1:100:end], yscale = :log10)
        # kmeaned = kmeans(reshape(log10.(residual), 1, :), 2)
        # θ3 = exp10(mean(kmeaned.centers))
        idx_blank = idx_blank[residual .< θ3]
        subsystem[idx_blank] .= id_subsys
        sets = sets[rand.(sets) .∉ Ref(idx_blank)]
        
        if sets |> isempty @info "sets exhausted"; break end
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


# function buck_recovery()
#     schedules = CSV.read("G:/DDM/bifurcation/buck_schedules.csv", DataFrame)[1:1:end,:]
#     vrbl = [:dV, :dI], [:V, :I]
#     cnfg = (; N = 1)
#     dt = 1e-7

#     @threads for dr in eachrow(schedules) # dr = eachrow(schedules)[1]
#         filename = "G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#         isfile(filename) && continue
#         data = CSV.read(replace(filename, "buck_rcvd" => "buck"), DataFrame);

#         add_subsystem!(data, vrbl, cnfg; θ1 = 1e+1, θ2 = 1e+0, min_rank = 2); # 30 sec
#         f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
#         Dtree = dryad(data, vrbl)

#         ic = collect(data[1, last(vrbl)])
#         try
#             ŷ = DataFrame(solve(f_, ic, dt, data.Vr, Dtree, data.Vr), last(vrbl))
#             CSV.write(filename, ŷ)
#         catch
#             @error "Error: $(lpad(dr.idx, 5, '0'))"
#         end
#     end
# end
# buck_recovery()

# _Vr = CSV.read("G:/DDM/bifurcation/buck/00001.csv", DataFrame).Vr
# idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
# @showprogress for dr in eachrow(schedules)
#     filename = "G:/DDM/bifurcation/buck_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#     isfile(filename) || continue
#     data = CSV.read(filename, DataFrame)

#     idx_sampled = diff(_Vr) .< 0
#     sampledV = data[Not(1), :V][idx_sampled]
#     append!(idcs, fill(dr.idx, length(sampledV)))
#     append!(hrzn, fill(dr.E, length(sampledV)))
#     append!(vrtc, sampledV)
# end
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1)
# CSV.write("G:/DDM/bifurcation/buck_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# png("G:/DDM/bifurcation/buck_recovered.png")


function soft_recovery()
schedules = CSV.read("G:/DDM/bifurcation/soft_schedules.csv", DataFrame)
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi, sign], λ = 1e-4)
cnfg = (; M = 3, λ = 1e-2)
dt = 1e-5

for dr in eachrow(schedules)[1:10:end,:] # dr = eachrow(schedules)[2001]
    filename = "G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) && continue
    data = CSV.read(replace(filename, "soft_rcvd" => "soft"), DataFrame);

    add_subsystem!(data, vrbl, cnfg; θ1 = 1e-8, θ2 = 1e-12, θ3 = 1e-5, min_rank = 21); # 15~
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
    Dtree = dryad(data, vrbl)
    
    try
        ic = collect(data[1, last(vrbl)])
        ŷ = DataFrame(solve(f_, ic, dt, data.t, Dtree), last(vrbl))
        CSV.write(filename, ŷ)
    catch
        @error "Error: $(lpad(dr.idx, 5, '0'))"
    end
end
end
soft_recovery()

idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr in eachrow(schedules)
    filename = "G:/DDM/bifurcation/soft_rcvd/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) || continue
    data = CSV.read(filename, DataFrame)

    idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
    sampledv = data[Not(1), :v][idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledv)))
    append!(hrzn, fill(dr.d, length(sampledv)))
    append!(vrtc, sampledv)
end
# CSV.write("G:/DDM/bifurcation/soft_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = .1);
# png("G:/DDM/bifurcation/soft_recovered.png")
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = .1, ylims = [-1, 1]);
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [430, 440])
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [440, 450], ticks = 440:450)
# scatter(idcs, vrtc, ms = 1, legend = :none, msw = 0, ma = 1, xlims = [550, 570])
# temp = DataFrame(; idcs, vrtc, hrzn)
# temp = temp[temp.idcs .∉ Ref([439, 448, 560]), :]
# scatter(temp.hrzn, temp.vrtc, ms = 1, legend = :none, msw = 0, ma = 1);
# png("temp")
# @info ""



# function hrnm_recovery()
#     schedules = CSV.read("G:/DDM/bifurcation/hrnm_schedules.csv", DataFrame)[1:1:end,:]
#     vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
#     cnfg = (; N = 3, f_ = [cos])
#     dt = 1e-3

#     for dr in eachrow(schedules) # dr = eachrow(schedules)[171]
#         filename = "G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#         isfile(filename) && continue
#         data = CSV.read(replace(filename, "hrnm_rcvd" => "hrnm"), DataFrame);

#         add_subsystem!(data, vrbl, cnfg; θ1 = 3e-2, θ2 = 1e-10, min_rank = 20); # 30 sec
#         f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
#         Dtree = dryad(data, vrbl)

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
#     filename = "G:/DDM/bifurcation/hrnm_rcvd/$(lpad(dr.idx, 5, '0')).csv"
#     isfile(filename) || continue
#     data = CSV.read(filename, DataFrame)

#     idx_sampled = abs.(diff(diff(data.z) ./ 1e-3)) .> 0.1
#     sampledx = data[Not(1, end), :x][idx_sampled]
#     append!(idcs, fill(dr.idx, length(sampledx)))
#     append!(hrzn, fill(dr.f, length(sampledx)))
#     append!(vrtc, sampledx)
# end
# scatter(hrzn, vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1)
# CSV.write("G:/DDM/bifurcation/hrnm_recovered.csv", DataFrame(; idcs, vrtc, hrzn))
# png("G:/DDM/bifurcation/hrnm_recovered.png")