include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

pool = CSV.read("G:/DDM/bifurcation/buck/02501.csv", DataFrame)
Random.seed!(0)
schedules = [DataFrame(idx = 1:100) DataFrame(rand(eachrow(pool), 100))]
for dr = eachrow(schedules)
    data = factory_buck(DataFrame, dr.idx, 40, ic = collect(dr)[[3,4]])
    CSV.write("G:/DDM/forecast/buck/$(lpad(dr.idx, 5, '0')).csv", data)
end
vrbl = [:dV, :dI], [:V, :I]
cnfg = (; N = 1)
dt = 1e-7
θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;
for dr in eachrow(schedules) # dr = eachrow(schedules)[end]
    filename = "G:/DDM/forecast/buck_prdt/$(lpad(dr.idx, 5, '0')).csv"
    isfile(filename) && continue
    data = CSV.read(replace(filename, "buck_prdt" => "buck"), DataFrame);

    add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 30 sec
    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)];
    Dtree = dryad(data, [:V, :I, :Vr])

    ic = collect(data[1, last(vrbl)])
    try
        ŷ = DataFrame(solve(f_, ic, dt, data.Vr, Dtree, data.Vr), last(vrbl))
        CSV.write(filename, ŷ)
    catch
        @error "Error: $(lpad(dr.idx, 5, '0'))"
    end
end


pool = CSV.read("G:/DDM/bifurcation/hrnm/00001.csv", DataFrame)
Random.seed!(0)
schedules = [DataFrame(idx = 1:100) DataFrame(rand(eachrow(pool), 100))]
# dr = first(eachrow(schedules))
for dr = eachrow(schedules)
    data = factory_hrnm(DataFrame, dr.idx, 0, ic = [0; collect(dr)[3:5]])
    CSV.write("G:/DDM/forecast/hrnm/$(lpad(dr.idx, 5, '0')).csv", data)
end
