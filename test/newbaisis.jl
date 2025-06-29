include("../core/header.jl")

import Base.rand
rand(df::AbstractDataFrame; n = 1) = df[rand(1:nrow(df), n), :]

data = factory_lorenz(DataFrame, 28, tspan = [0, 100])[1:100:end, :]
vrbl = [:dx, :dy, :dz], [:x, :y, :z]
strv = string.(last(vrbl))
XY = rand(data, n = 100)
X = Matrix(XY[:, last(vrbl)])
Y = Matrix(XY[:, first(vrbl)])


cnfg = cook(last(vrbl); poly = [0, 1, 3])
# cnfg = cnfg[Not([4,6,24,1]), :]
ΘX = Θ(X, cnfg)
f = SINDy(XY, vrbl, cnfg, λ = 1e-2)
print(f)
Θ(X, cnfg)
