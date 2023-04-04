@time using DataDrivenDiffEq
@time using DataDrivenSparse

@time using CSV, DataFrames, Dates
function fillmissing!(data)
    for col in eachcol(data)
        while true
            bit_missing = ismissing.(col)
            if sum(bit_missing) == 0 break end
            col[bit_missing] .= ((circshift(col, -1) + circshift(col, 1)) / 2)[bit_missing]
            bit_missing = ismissing.(col)
            col[bit_missing] .= circshift(col, 1)[bit_missing]
        end
    end
end

data = DataFrame()
DATA_ = []
for fn = readdir("data")
    itemname = Symbol(first(split(fn, " ")))
    if (itemname == :알루미늄) || (itemname == :오렌지) continue end
    push!(DATA_, CSV.read("data/" * fn, DataFrame))
    DATA_[end].날짜 = Date.(DATA_[end].날짜, dateformat"y- m- d")
    select!(DATA_[end], ["날짜", "종가"])
    rename!(DATA_[end], [:t, itemname])
    if eltype(DATA_[end][:, 2]) <: AbstractString
        DATA_[end][:, 2] = replace.(DATA_[end][:, 2], "," => "")
        DATA_[end][!, 2] = parse.(Float64, DATA_[end][:, 2])
    end
    if isempty(data)
        data = deepcopy(DATA_[end])
    else
        leftjoin!(data, DATA_[end], on = :t, makeunique = true)
    end
end
data = data[data.t .< Date(2022, 1, 1), :]
sort!(data, :t)
fillmissing!(data)
@assert data == dropmissing(data) # 데이터 무결성 검사
dropmissing!(data)
trng = data[data.t .< Date(2021, 1, 1), :]
test = data[data.t .≥ Date(2020, 12, 31), :]

selected = [:WTI유, :구리, :금]

# scatter(trng[:, :WTI유], trng[:, :천연가스])


# X = Matrix(data[:, Not(:t)])'
X = Matrix(trng[:, selected])'
DX = diff(X, dims = 2); X = X[:, 1:(end-1)]
ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)

opt = STLSQ(10^(-4), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq)
P = get_parameter_values(soleq)

@time using OrdinaryDiffEq
u0 = X[:, end]
tspan = (30, 50)
dt = 1
DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)

@time using Plots
p1 = plot(DDM, idxs = (0,1))
p2 = plot(DDM, idxs = (0,2))
p3 = plot(DDM, idxs = (0,3))
plot!(p1, 0:30, X[1, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "WTI")
plot!(p2, 0:30, X[2, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "Copper")
plot!(p3, 0:30, X[3, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "Gold")
plot!(p1, 30:50, test[1:21, 2], color = :black, label = "Data")
plot!(p2, 30:50, test[1:21, 3], color = :black, label = "Data")
plot!(p3, 30:50, test[1:21, 4], color = :black, label = "Data")
plot(p1, p2, p3, layout = (3, 1), size = (800, 600))
