include("../src/factorio.jl")

d_range = 0.1:0.0001:0.3
schedule = DataFrame(idx = eachindex(d_range), d = d_range)

## -------------- begin bifurcation diagram -------------- ##
# xdots = Float64[]; ydots = Float64[]
# for dr in ProgressBar(eachrow(schedule))
#     data = factory_soft(dr.idx, dr.d)
    
#     idx = [false; diff(abs.(data.u) .> (dr.d/2)) .> 0]

#     sampled = data.v[idx]
#     append!(xdots, fill(dr.d, length(sampled)))
#     append!(ydots, sampled)
# end
# @time a1 = scatter(xdots, ydots,
# xlabel = L"d", ylabel = "Impact velocity",
# label = :none, msw = 0, color = :black, ms = 0.5, alpha = 0.5, size = (700, 300))
# png(a1, "soft_bifurcation")
## -------------- end bifurcation diagram -------------- ##

include("../src/ML.jl")

using Plots, LaTeXStrings; default(msw = 0)

dr = eachrow(schedule)[501]
data = factory_soft(dr.idx, dr.d)
_data = data[1:1:end, :]
points = col_normalize([_data.u _data.v _data.dv])'
a2 = scatter(eachrow(points)...);

using Clustering
ε = 0.0001
@time dbscaned = dbscan(points, ε)
dbscaned.counts
subsystem = dbscaned.assignments
a3 = scatter(eachrow(points)..., title = L"ε=" * "$ε", color = subsystem, label = :none, ms = 1);
png(a3, "230905 2")

__data = deepcopy(_data)
__data = _data[:, Not([:dt, :du])]
__data[:, :label] = subsystem
CSV.write("data/softimpact15.csv", __data)

include("../src/DDM.jl")

gdf_ = groupby(__data, :label)
mtrx_sindy = []
for gdf in gdf_
    X = poly_basis(col_func(Matrix(gdf[:, 1:3]), [cospi, sinpi, abs, sign]), 2)
    Y = Matrix(gdf[:, 3:4])
    push!(mtrx_sindy, STLSQ(X, Y, 0.1))
end
mtrx_sindy[1]
mtrx_sindy[2]
mtrx_sindy[3]
mtrx_sindy[4]

mtrx_sindy[1] - mtrx_sindy[3]

gdf = gdf_[1]