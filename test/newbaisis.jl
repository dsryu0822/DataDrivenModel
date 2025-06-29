include("../core/header.jl")

import Base.rand
rand(df::AbstractDataFrame; n = 1) = df[rand(1:nrow(df), n), :]

data = factory_lorenz(DataFrame, 28, tspan = [0, 100])[1:100:end, :]
vrbl = [:dx, :dy, :dz], [:x, :y, :z]
strv = string.(last(vrbl))
XY = rand(data, n = 100)
X = Matrix(XY[:, last(vrbl)])
Y = Matrix(XY[:, first(vrbl)])
N = 3
function factorsort(xyz)
    dict_sup = Dict(1 => "", 2 => "²", 3 => "³", 4 => "⁴", 5 => "⁵", 6 => "⁶", 7 => "⁷", 8 => "⁸", 9 => "⁹")
    xyz_ = split(xyz, '⋅')
    xyz1 = sort(unique(xyz_))
    deg_ = [count(xyz_ .== x) for x in xyz1]
    return join(xyz1 .* [dict_sup[d] for d in deg_], '⋅')
end


constant(x) = 1
recipe = DataFrame(index = Int64[], deg = Int64[], term = String[], func = [], funh = [], vecv = [])
idx = 0
# push!(recipe, [1, 0, "1", constant, [1]])
for i in 1:size(X, 2)
    idx += 1
    push!(recipe, [idx, 1, strv[i], identity, vec, [i]])
end
for d = 1:N
    for i in eachindex(strv)
        rcp = recipe[recipe.deg .== d, :]
        for j in i:nrow(rcp)
            idx += 1
            push!(recipe, [idx, d + 1, "$(strv[i])⋅$(rcp.term[j])", prod, eachrow, sort([i; rcp.vecv[j]])])
        end
    end
end
pushfirst!(recipe, [0, 0, "1", constant, vec, [1]])
recipe.term = factorsort.(recipe.term)
push!(recipe, [67, 0, "cos x", cos, vec, [1]])

dr = eachrow(recipe)[19]

dr.func.(eachrow(X[:, dr.vecv]))
dr.func.(vec(X[:, dr.vecv]))

Θ(X, recipe) = hcat([dr.func.(dr.funh(X[:, dr.vecv])) for dr in eachrow(recipe)]...)

ΘX = Θ(X, recipe)
ΘX \ Y