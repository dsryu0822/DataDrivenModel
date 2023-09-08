include("../src/factorio.jl")
include("../src/ML.jl")
include("../src/DDM.jl")


using Plots, LaTeXStrings; default(msw = 0, color = :black)

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx = eachindex(d_range), d = d_range)

cd("//155.230.155.221/ty/DS")

dr = eachrow(plan)[501]
data = factory_soft(dr.idx, dr.d)

data[:, :ddv] = [diff(data.dv)[1]; diff(data.dv)]
_data = data[1:1:end, :]

# histogram(diff(__data.dv), yscale = :log10, xlabel = "u'''"); png("230907 2");
# scatter(__data.u[Not(1)], diff(__data.dv), xlabel = "u", ylabel = "u'''"); png("230907 4");

# histogram(log.(abs.(_data.ddv)), xlabel = "log.(abs.(_data.ddv))", yscale = :log10); png("230907 2")

# points = col_normalize([_data.u _data.v _data.dv])'
histogram(log.(abs.(_data.ddv)), xlabel = "log.(abs.(_data.ddv))", yscale = :log10); png("230907 3")

subsystem = (log.(abs.(_data.ddv)) .< -10) .+ 1
scatter(_data.u, _data.v, _data.dv, color = subsystem, label = :none, ms = 1, alpha = 0.5); png("230907")

valnames = ["t" "u" "v" "cos(t)" "cos(u)" "cos(v)" "abs(t)" "abs(u)" "abs(v)" "sign(t)" "sign(u)" "sign(v)"] |> vec
# valnames = poly_basis(valnames', 2, forcing = true) |> vec

_data[:, :subsystem]  = subsystem

gdf_ = groupby(_data, :subsystem)
mtrx_ = []
for gdf in gdf_
    X = poly_basis(Matrix(gdf[:, [:t, :u, :v]]), 3)
    # X = col_func(Matrix(gdf[:, [:t, :u, :v]]), [cospi, abs, sign])
    Y = Matrix(gdf[:, [:du, :dv]])
    push!(mtrx_, STLSQ(X, Y, 0.1))
end
mtrx_[1]
mtrx_[2]
[valnames mtrx_[1]]
[valnames mtrx_[2]]

scatter(cospi.(gdf_[1].t), gdf_[1].dv)

Plots.gr()
scatter(gdf_[1].u, gdf_[1].v, gdf_[1].dv, zlims = (-150, 200), ms = 0.1); png("temp1")
scatter(gdf_[2].u, gdf_[2].v, gdf_[2].dv, zlims = (-150, 200), ms = 0.1); png("temp2")


STLSQ(
    [cospi.(gdf_[1].t) gdf_[1].u gdf_[1].v sign.(gdf_[1].u)],
    [gdf_[1].du gdf_[1].dv], 0.1)

z = [cospi.(gdf_[1].t) gdf_[1].u gdf_[1].v sign.(gdf_[1].u)]

points = reshape(log.(abs.(_data.ddv)), 1, :)
ε = 0.1
@time dbscaned = dbscan(points, ε); dbscaned.counts
subsystem = dbscaned.assignments
# scatter(_data.u, _data.v, _data.dv, title = L"ε=" * "$ε", color = subsystem, label = :none, ms = 1);