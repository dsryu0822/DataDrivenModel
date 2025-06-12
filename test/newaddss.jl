include("../core/header.jl")

using JLD2
using Graphs

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5; tspan = [0, 20]; #θ = 1e-16;
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
    test = CSV.read("hyperparameter/0 test.csv", DataFrame)
else
    # data = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    # test = factory_soft(DataFrame, 0.1, ic = [0, .000129715, 0.301469]; tspan, dt)
    temp = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    data = temp[1:(nrow(temp) ÷ 2), :]
    test = temp[(nrow(temp) ÷ 2):end, :]
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
    CSV.write("hyperparameter/0 test.csv", test, bom = true)
end

import Base.rand
rand(df::AbstractDataFrame; n = 1) = df[rand(1:nrow(df), n), :]
doublerange(n) = Base.product(1:n, 1:n)
centerpick(n, m) = round.(Int64, range(1, n, m+2))[2:end-1]
has_neighbor(arr) = !isempty(intersect(arr, [arr .+ 1; arr .- 1]))

jumpt = detect_jump(data, vrbl)
sets = set_divider(jumpt)
datasets = [data[set, :] for set in sets]
ground_truth = [ifelse(rand(ds.u) > 0.05, 3, ifelse(rand(ds.u) < -0.05, 2, 1)) for ds in datasets]

N_ = 0:3
M_ = 0:3
schedules = DataFrame(ID = Int64[], N = Int64[], M = Int64[], worked = String[], ss = Int64[], totalerror = Float64[])
for (N, M) in Base.product(N_, M_) push!(schedules, (0, N, M, "", 0, 0.0, )) end
schedules.ID = 1:size(schedules, 1)

for dr = eachrow(schedules)
    ID = dr.ID; N = dr.N; M = dr.M
    cnfg = (; N, M, λ = 0)

    # cnfg = (; N = 1, M = 1, λ = 0)
    # # min_m = minimum(nrow.(datasets))
    n = length(datasets)
    min_m = SINDy(datasets[1][1:10, :], vrbl...; cnfg...).sparse_matrix.m
    # sampled = rand.(datasets; n = min_m)
    sampled = [ds[centerpick(nrow(ds), min_m), :] for ds in datasets]
    f__ = fill(SINDy(datasets[1][1:10, :], vrbl...; cnfg...), 29, 29)
    mse__ = fill(Inf, n, n)
    dist__ = zeros(n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end
    f_ = f__[argmin(mse__, dims = 2)]
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        dist__[i, j] = sum(abs2, f_[i].sparse_matrix - f_[j].sparse_matrix)
        dist__[j, i] = dist__[i, j]
    end

    nzθ = vec(filter(!iszero, dist__))
    candy = DataFrame(
        delta = nzθ,
        ss = [length(connected_components(SimpleGraph(dist__ .< θ))) for θ in nzθ]
    )
    sort!(candy, :delta, rev = true)
    unique!(candy, :ss, keep = :last)
    θ_ = candy.delta
    candy.delta .= sqrt.([10candy.delta[1]; candy.delta[1:end-1]] .* candy.delta)

    sargs = (; legend = :none, msw = 0, ms = 2)
    sdelta = scatter(title = "N = $N, M = $M", ylabel = L"\delta_{ij}", size = [600, 1200], left_margin = 10mm)
    scatter!(sdelta, nzθ, yscale = :log10; sargs...)
    hline!(sdelta, candy.delta, color = :gray)
    scatter!(sdelta, fill(length(nzθ) + 50, nrow(candy)), candy.delta, text = string.(candy.ss), ms = 0)
    png(sdelta, "hyperparameter/$ID d n=$(N)_m=$(M).png")

    num_subsys = []
    totalerror = []
    bit_feasible = []
    for θ in θ_
        subsys = connected_components(SimpleGraph(dist__ .< θ))
        push!(num_subsys, length(subsys))
        push!(bit_feasible, !(any(length.(subsys) .== 1) || any(has_neighbor.(subsys))))
        
        label = zeros(Int64, length(sampled))
        for s in eachindex(subsys) label[subsys[s]] .= s end

        ref_ = [SINDy([sampled[ss]...;], vrbl...; cnfg...) for ss in subsys]
        erroref_ = [[] for _ in eachindex(ref_)]
        for k in eachindex(f_)
            push!(erroref_[label[k]], norm(f_[k].sparse_matrix - ref_[label[k]].sparse_matrix))
        end
        push!(totalerror, mean(mean.(erroref_)))
    end
    argθ = argmin(totalerror + (Inf .* .!bit_feasible))
    θ = θ_[argθ]
    dr .= [ID, N, M, ifelse(any(bit_feasible), "✅", "❌"), num_subsys[argθ], totalerror[argθ]]
    # scatter(num_subsys, θ_, yscale = :log10, xticks = 1:length(num_subsys))
    scatter(num_subsys, totalerror, yscale = :log10, xticks = 1:length(num_subsys),
            shape = ifelse.(bit_feasible, :o, :x), legend = :none,
            ylabel = "Total equation error", xlabel = "Number of subsystems", title = "N = $N, M = $M",)
    png("hyperparameter/$ID n=$(N)_m=$(M).png")
    open("hyperparameter/$ID n=$(N)_m=$(M).txt", "w") do scroll
        A = (dist__ .< θ)
        subsys = connected_components(SimpleGraph(A))
        f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg..., λ = 1e-3) for ss in subsys]
        println(scroll, "ID = $ID, N = $N, M = $M")
        println(scroll, "Number of subsystems: $(length(subsys))")
        for (k, f) in enumerate(f_)
            print(scroll, pretty_table(f))
            println(scroll, "MSE = $(f.MSE)\n")
        end
        println(scroll, "Total equation error: $(dr.totalerror)")
    end
end
schedules

lperform = plot()
for df in groupby(schedules, :M)
    plot!(lperform, df.N, df.totalerror, label = "M = $(df.M[1])", lw = 2)
end
plot(lperform, xlabel = "N", yscale = :log10, legend = :topleft, size = [400, 400], ylabel = "totral equation error E")


for k in eachindex(f_)
    open("scroll.txt", "w") do scroll
        println(scroll, "D$k, number of datapoint: $(nrow(datasets[k]) + 3)")
        print(scroll, pretty_table(f_[k]))
        println(scroll, "MSE = $(f_[k].MSE)\n")
    end
end

function labeling4!(data, vrbl, cnfg; dos = 0, L = 0)
    jumpt = detect_jump(data, vrbl)

    sets = set_divider(jumpt)
    datasets = [data[set, :] for set in sets]
    # ground_truth = [ifelse(rand(ds.u) > 0.05, 3, ifelse(rand(ds.u) < -0.05, 2, 1)) for ds in datasets]

    if L |> iszero
        L = zeros(length(sets), length(sets))
        for i in eachindex(sets), j in eachindex(sets)
            if i < j
                # print("($i, $j)")
                data_alien = [datasets[[i, j]]...;]
                f = SINDy(data_alien, vrbl...; cnfg..., λ = 0)
                L[i, j] = f.MSE
            end
        end
        L = Symmetric(L)
    end

    label = zeros(Int64, length(sets))
    idcs = eachindex(sets)
    f_ = []
    # tol_ = []
    subsys = 0
    mse_ = []
    for subsys in eachindex(sets) # subsys += 1
        idx_longest = argmax(iszero.(label) .* length.(sets))
        idx_farthest = idcs[argmax(L[idx_longest, :][idcs])]
        bit_inside = L[idx_longest, :] .< L[idx_farthest, :]
        candidates = findall(bit_inside) ∩ idcs
        inside = candidates[L[idx_longest, candidates] .< minimum(L[idx_longest, setdiff(idcs, candidates)])]
        outside = setdiff(idcs, inside)

        if length(outside) ≤ 2 ||
           maximum(L[idcs, idcs]) < maximum([mse_; 0])
            label[iszero.(label)] .= subsys
        else
            label[inside] .= subsys
            # push!(tol_, minimum(L[idx_longest, outside]))
        end
        push!(f_, SINDy([datasets[label .== subsys]...;], vrbl...; cnfg..., λ = 0))
        push!(mse_, f_[end].MSE)
        
        idcs = setdiff(idcs, findall(label .== subsys))
        if isempty(idcs) break end
    end
    

        #   idx_longest = argmax(iszero.(label) .* length.(sets))
#         idx_farthest = idcs[argmax(L[idx_longest, :][idcs])]
#         bit_inside = L[idx_longest, :] .< L[idx_farthest, :]
#         inside = findall(bit_inside) ∩ idcs
#         outside = setdiff(idcs, inside)
# # SINDy([datasets[iszero.(label)]...;], vrbl...; cnfg..., λ = 0).MSE
#         candidates = inside[L[idx_longest, inside] .< minimum(L[idx_longest, outside])]

#         if length(idcs) == (length(candidates) + 1) ||
#             maximum(L[idcs, idcs]) < maximum([mse_; 0])
#             label[iszero.(label)] .= subsys
#         else
#             label[candidates] .= subsys
#         end
#         push!(f_, SINDy([datasets[label .== subsys]...;], vrbl...; cnfg..., λ = 0))
#         # push!(mse_, f_[end].MSE)
        
#         idcs = setdiff(idcs, findall(label .== subsys))
#         if isempty(idcs) break end
    # idcs
    # [ground_truth label]

    data[!, :label] = zeros(Int64, nrow(data))
    for (k, lbl) in enumerate(label)
        data.label[sets[k]] .= lbl
    end
    bit_zero = iszero.(data.label)
    # f_ = [SINDy(data, vrbl...; cnfg...) for data in groupby(data[.!bit_zero,:], :label)]
    data.label[bit_zero] .= argmin.(eachrow(stack([sum.(abs2, f.(eachrow(data[bit_zero, last(vrbl)])) - collect.(eachrow(data[bit_zero, first(vrbl)]))) for f in f_])))

    return data
end
# labeling4!(data, vrbl, cnfg)

N_ = 0:3
M_ = 0:3
schedules = DataFrame(ID = Int64[], N = Int64[], M = Int64[], ss = Int64[], trngSSE = Float64[], testSSE = Float64[], runtime = [])
for (N, M) in Base.product(N_, M_) push!(schedules, (0, N, M, 0, 0.0, 0.0, 0)) end
schedules.ID = 1:size(schedules, 1)

if !isfile("hyperparameter/0 L_.jld2")
    L_ = []
    for dr = eachrow(schedules)
        ID = dr.ID; λ = 0; n = dr.N; m = dr.M
        cnfg = (; N = n, M = m, λ = λ)
        
        jumpt = detect_jump(data, vrbl; dos = 0)

        sets = set_divider(jumpt)
        datasets = [data[set, :] for set in sets]
        # ground_truth = [ifelse(rand(ds.u) > 0.05, 3, ifelse(rand(ds.u) < -0.05, 2, 1)) for ds in datasets]

        L = zeros(length(sets), length(sets))
        for i in eachindex(sets), j in eachindex(sets)
            if i < j
                # print("($i, $j)")
                data_alien = [datasets[[i, j]]...;]
                f = SINDy(data_alien, vrbl...; cnfg..., λ = 0)
                L[i, j] = f.MSE
            end
        end
        L = Symmetric(L)
        push!(L_, L)
        JLD2.save("hyperparameter/0 L_.jld2", "L_", L_)
    end
else
    L_ = JLD2.load("hyperparameter/0 L_.jld2")["L_"]
end

for dr = eachrow(schedules)
    tic = now()
    try
        ID = dr.ID; λ = 0; n = dr.N; m = dr.M
        cnfg = (; N = n, M = m, λ = λ)
        
        trng = deepcopy(data)
        labeling4!(trng, vrbl, cnfg; L = L_[ID])
        if 0 ∈ trng.label
            dr.ss = -1
        else
            dr.ss = length(unique(trng.label))
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :label)] # print.(f_)
            Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
            pred1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, trng.t, Dtree), last(vrbl))
            dr.trngSSE = sum(abs2, trng.u - pred1.u)
            pred2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, test.t, Dtree), last(vrbl))
            dr.testSSE = sum(abs2, test.u - pred2.u)

            open("hyperparameter/17 equation.txt", "a") do io
                println(io, "ID = $(dr.ID), N = $(dr.N), M = $(dr.M)")
                for f in f_
                    print(io, pretty_table(f))
                end
                println(io, "---------------------------------")
            end
        
            CSV.write("hyperparameter/$ID n=$(n)_m=$(m).csv", pred2, bom = true)
        end
    catch
        dr.trngSSE = -1
        dr.testSSE = -1
    end
    dr.runtime = (now() - tic).value / 60000
    println(schedules)
    CSV.write("hyperparameter/17 result.csv", schedules, bom = true)
end

data[:, :label] = [ifelse(dr.u > 0.05, 3, ifelse(dr.u < -0.05, 2, 1)) for dr in eachrow(data)]
count(data.label .== 1) / nrow(data) # 0.5
count(data.label .== 2) / nrow(data) # 0.5
count(data.label .== 3) / nrow(data) # 0.5
plt = plot(data.t, data.u, color = data.label, label = "u", lw = 2, dpi = 300, size = [1920, 1080]);
png(plt, "temp.png")

[["D$n" for n in 1:29] nrow.(datasets) .+ 3]
cnfg = (; N = 1, M = 1, λ = 1e-1)
f_ = [SINDy(ds, vrbl...; cnfg...) for ds in datasets]
for k in 1:length(f_)
    println("\nD$k, number of datapoint: $(nrow(datasets[k])+3)")
    print(f_[k])
    println("MSE = $(f_[k].MSE)")
end
rank(Θ(datasets[1]; cnfg...))

[1, 3, 5, 9, 11, 12, 13, 15, 19, 21, 22, 25, 28, 29]
1:2:4 |> typeof
lasticds = cumsum(nrow.([d[k, :] for (d,k) in zip(datasets, StepRange.(1, nrow.(datasets) .÷ 100, nrow.(datasets)))]))
expanded = [[d[k, :] for (d,k) in zip(datasets, StepRange.(1, nrow.(datasets) .÷ 100, nrow.(datasets)))]...;]
plt = scatter(expanded.u, color = expanded.label, xticks = (lasticds, 1:29), msw = 0, ms = 1, dpi = 300, legend = :none, size = 50 .* [16, 9]);
png(plt, "temp.png")



f_ = [SINDy(ds, vrbl...; cnfg...) for ds in [d[k, :] for (d,k) in zip(datasets, StepRange.(1, nrow.(datasets) .÷ 100, nrow.(datasets)))]]
for k in 1:length(f_)
    println("\nD$k, number of datapoint: $(nrow(datasets[k])+3)")
    print(f_[k])
    println("MSE = $(f_[k].MSE)")
end


scatter(vec(mse__), yscale = :log10, title = "N = 3, M = 3", ylabel = "MSE: " * L"S_i \cup S_j", msw = 0, ms = 2, legend = :none, ylims = [-Inf, 0], yticks = exp10.(-30:10:0))
connected_components(SimpleGraph(mse__ .< 1e-25))
maximum(mse__[mse__ .< Inf])