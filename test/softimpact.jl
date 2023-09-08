include("../src/factorio.jl")
include("../src/ML.jl")
include("../src/DDM.jl")


using Plots, LaTeXStrings; default(msw = 0, color = :black)

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx = eachindex(d_range), d = d_range)

## -------------- begin bifurcation diagram -------------- ##
# xdots = Float64[]; ydots = Float64[]
# for dr in ProgressBar(eachrow(plan))
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


# dr = eachrow(plan)[501]
# data = factory_soft(dr.idx, dr.d)
# _data = data[1:1:end, :]
# points = col_normalize([_data.u _data.v _data.dv])'
# a2 = scatter(eachrow(points)...);

# using Clustering
# ε = 0.0001
# @time dbscaned = dbscan(points, ε); dbscaned.counts
# subsystem = dbscaned.assignments
# a3 = scatter(eachrow(points)..., title = L"ε=" * "$ε", color = subsystem, label = :none, ms = 1);
# png(a3, "230905 2")

# __data = deepcopy(_data)
# __data = _data[:, Not([:dt, :du])]
# __data[:, :label] = subsystem
# CSV.write("data/softimpact15.csv", __data)


# __data = CSV.read("data/softimpact15.csv", DataFrame)

# gdf_ = groupby(__data, :label) # gdf = gdf_[1]
# mtrx_sindy = []
# for gdf in gdf_
#     X = poly_basis(col_func(Matrix(gdf[:, 1:3]), [cospi, abs, sign]), 2)
#     Y = Matrix(gdf[:, 3:4])
#     push!(mtrx_sindy, STLSQ(X, Y, 0.01))
# end
# mtrx_sindy[1]
# mtrx_sindy[2]
# mtrx_sindy[3]
# mtrx_sindy[4]

# valnames = poly_basis(["t" "u" "v" "cost" "cosu" "cosv" "abst" "absu" "absv" "signt" "signu" "signv"], 2, forcing = true) |> vec

# [valnames mtrx_sindy[3]]
# 0.112473 + 0.443764 + 0.443764

# X = poly_basis(col_func(Matrix(gdf[:, 1:3]), [cospi, abs, sign]), 2)
# det(X'X)

# gdf = gdf_[1][1:10:end, :]
# STLSQ([cospi.(gdf.t) gdf.v], Matrix(gdf[:, 3:4]), 0.01)
# scatter(cospi.(gdf.t), gdf.dv, color = 1); png("230907 2")

# using Clustering
# ε = 0.00005
# @time dbscaned = dbscan(points, ε, min_cluster_size = 91); dbscaned.counts
# subsystem = dbscaned.assignments
# a3 = scatter(eachrow(points)..., title = L"ε=" * "$ε", color = subsystem, label = :none, ms = 1);
# png(a3, "230906 1")

cd("//155.230.155.221/ty/DS")

dr = eachrow(plan)[501]
data = factory_soft(dr.idx, dr.d)

data[:, :ddv] = [eps(); diff(data.dv)]
_data = data[1:1:end, :]

# histogram(diff(__data.dv), yscale = :log10, xlabel = "u'''"); png("230907 2");
# scatter(__data.u[Not(1)], diff(__data.dv), xlabel = "u", ylabel = "u'''"); png("230907 4");

# histogram(log.(abs.(_data.ddv)), xlabel = "log.(abs.(_data.ddv))", yscale = :log10); png("230907 2")

# points = col_normalize([_data.u _data.v _data.dv])'
points = reshape(log.(abs.(_data.ddv)), 1, :)
ε = 0.1
@time dbscaned = dbscan(points, ε); dbscaned.counts
subsystem = dbscaned.assignments



scatter(_data.u, _data.v, _data.dv, title = L"ε=" * "$ε", color = subsystem, label = :none, ms = 1);
png(a3, "230907 3")

histogram(vec(points))
scatter(_data.u, _data.v, _data.dv, color = (abs.(_data.ddv) .> 5) .+ 1, label = :none, ms = 1); png("temp")
