include("../core/header.jl")
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-4; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;

data = CSV.read("lyapunov/hrnm_traj/00001.csv", DataFrame)

f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)

typeof(s)
function functionalizer(s::STLSQresult)
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.matrix, dims = 1))
    return v -> substitute(fnexp, Dict(rname .=> v))
end
g_ = functionalizer.(f_)

using BenchmarkTools

x = collect(data[5, 1:4])
@btime f_[1]($x)
@btime g_[1]($x)


f_