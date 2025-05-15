include("../core/header.jl")

using JLD2

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5; tspan = [0, 10]; #θ = 1e-16;
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
    test = CSV.read("hyperparameter/0 test.csv", DataFrame)
else
    data = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    test = factory_soft(DataFrame, 0.1, ic = [0, .000129715, 0.301469]; tspan, dt)
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
    CSV.write("hyperparameter/0 test.csv", test, bom = true)
end

cnfg = (; N = 3, M = 1, λ = 0)

function labeling4!(data, vrbl, cnfg; dos = 0, L = 0)
jumpt = detect_jump(data, vrbl; dos = 0)

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

for subsys in eachindex(sets) # subsys = 3
        # if label |> !iszero
        #     # worst_mse = SINDy([datasets[setdiff(eachindex(sets), idcs)]...;], vrbl...; cnfg..., λ = 0).MSE
        #     # if maximum(L[idcs, idcs]) < worst_mse
        #     # if maximum(L[idcs, idcs]) < minimum(maximum(L[idcs, setdiff(eachindex(sets), idcs)], dims = 2))
        #     if maximum(L[idcs, idcs]) < minimum(L[idcs, setdiff(eachindex(sets), idcs)])
        #         label[idcs] .= subsys
        #         break
        #     end
        # end

        idx_longest = argmax(iszero.(label) .* length.(sets))
        idx_farthest = idcs[argmax(L[idx_longest, :][idcs])]
        bit_inside = L[idx_longest, :] .< L[idx_farthest, :]
        inside = findall(bit_inside) ∩ idcs
        outside = setdiff(idcs, inside)

        candidates = inside[L[idx_longest, inside] .< minimum(L[idx_longest, outside])]
        if length(idcs) == (length(candidates) + 1)
            label[iszero.(label)] .= subsys
            break
        else
            label[candidates] .= subsys
            idcs = setdiff(idcs, candidates)
        end
    end
    # idcs
    # [ground_truth label]

    data[!, :label] = zeros(Int64, nrow(data))
    for (k, lbl) in enumerate(label)
        data.label[sets[k]] .= lbl
    end
    bit_zero = iszero.(data.label)
    f_ = [SINDy(data, vrbl...; cnfg...) for data in groupby(data[.!bit_zero,:], :label)]
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
            pred1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, first(tspan):dt:last(tspan), Dtree)[1:(end-1),:], last(vrbl))
            dr.trngSSE = sum(abs2, trng.u - pred1.u)
            pred2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, first(tspan):dt:last(tspan), Dtree)[1:(end-1),:], last(vrbl))
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
