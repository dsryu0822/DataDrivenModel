@time using DataDrivenDiffEq
@time using CSV, DataFrames, Dates
@time using DataDrivenSparse


data = DataFrame()
DATA_ = []
for fn = readdir("data")
    push!(DATA_, CSV.read("data/" * fn, DataFrame))
    DATA_[end].날짜 = Date.(DATA_[end].날짜, dateformat"y- m- d")
    select!(DATA_[end], ["날짜", "종가"])
    rename!(DATA_[end], [:t, Symbol(first(split(fn, " ")))])
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
dropmissing!(data)

X = Matrix(data[:, Not(:t)])'
DX = X - circshift(X, (0,1))

ddprob = ContinuousDataDrivenProblem(X, DX)

# @variables u1 u2 u3
# u = [u1; u2; u3]
@variables u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12
u = [u1; u2; u3; u4; u5; u6; u7; u8; u9; u10; u11; u12]
basis = Basis(polynomial_basis(u, 2), u)

opt = STLSQ(exp10.(-7:5))
# https://github.com/SciML/DataDrivenDiffEq.jl/blob/cae0c79e0d7b9eef48f9350e01f69b0798b447bb/lib/DataDrivenSparse/src/algorithms/STLSQ.jl#L97-L120
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol)
println(soleq) # soleq[1].rhs
is_converged(ddsol)

