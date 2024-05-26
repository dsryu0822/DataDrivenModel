include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

pool = CSV.read("G:/DDM/cached_buck.csv", DataFrame)
Random.seed!(0)
schedules = [DataFrame(idx = 1:100) DataFrame(rand(eachrow(pool), 100))]
for dr = eachrow(schedules)
    data = factory_buck(DataFrame, dr.idx, 40, ic = collect(dr)[[3,4]])
    CSV.write("G:/DDM/forecast/buck_test/$(lpad(dr.idx, 5, '0')).csv", data)
end
vrbl = [:dV, :dI], [:V, :I]
cnfg = (; N = 1)
dt = 1e-7
θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;

add_subsystem!(pool, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(pool, :subsystem)];
Dtree = dryad(pool, [:V, :I, :Vr])

for dr = eachrow(schedules)
    filename = "G:/DDM/forecast/buck_prdt/$(lpad(dr.idx, 5, '0')).csv"
    ic = collect(dr[last(vrbl)])
    try
        ŷ = DataFrame(solve(f_, ic, dt, data.Vr, Dtree, data.Vr), last(vrbl))
        CSV.write(filename, ŷ)
    catch
        @error "Error: $(lpad(dr.idx, 5, '0'))"
    end
end

buck_result = []
@showprogress for dr = eachrow(schedules)
    prdt = CSV.read("G:/DDM/forecast/buck_prdt/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    test = CSV.read("G:/DDM/forecast/buck_test/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    add_subsystem!(test, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...);
    idx_cng = findall(.!iszero.(diff(test.subsystem)))
    push!(buck_result, abs(test.V[idx_cng[10]] - prdt.V[idx_cng[10]]) < 1e-2)
end
count(buck_result)