using LinearAlgebra, Clustering, CSV, DataFrames

function col_normalize(M)
    return M ./ norm.(eachcol(M))'
end


# DATA = CSV.read("G:/buck/buck_000006.csv", DataFrame)[end-100000:end,:]
# Y = select(DATA, [:dV, :dI]) |> Matrix
# X = select(DATA, [ :V,  :I]) |> Matrix
# XY = [X Y]

# dbs = dbscan(col_normalize(Y)', 0.001); nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
# using Plots
# scatter(eachcol(Y)..., msw = 0, color = dbs.assignments, label = :none)

using Flux
struct Fourier0
    an
    bn
    L
end
Fourier0(m::Integer) = Fourier0(randn(Float32, m), randn(Float32, m), Float32[0.0])
function (layer::Fourier0)(x)
    m = length(layer.an)
    nxL⁻¹ = (1:m) * x[end, :]' .* exp(layer.L[1])
    return [x[1:end-1, :]; sum(
        [cospi.(nxL⁻¹) .* layer.an
       ; sinpi.(nxL⁻¹) .* layer.bn], dims = 1)]
end
Flux.@functor Fourier0

function Base.show(io::IO, l::Fourier0)
    print(io, "Fourier(", length(l.an), ")")
end
# ------------------------------

struct Fourier
    m::Int64
    an::Array{Float32, 1}
    bn::Array{Float32, 1}
    L::Float32
end
function Base.show(io::IO, l::Fourier)
    print(io, "Fourier(", l.m, ")")
end
Fourier(m::Integer) = Fourier(m, randn(Float32, m), randn(Float32, m), zero(Float32))

function (layer::Fourier)(x)
    nxL⁻¹ = (1:layer.m) .* last(eachrow(x))' * exp(layer.L)
    return [x[1:end-1, :]; sum(
        [cospi.(nxL⁻¹) .* layer.an
       ; sinpi.(nxL⁻¹) .* layer.bn], dims = 1)]

    # nxL⁻¹ = (1:layer.m) .* last(eachrow(x))' * exp(layer.L)
    # return [x[1:end-1, :]; sum([sinpi.(nxL⁻¹); cospi.(nxL⁻¹)] .* layer.coeff, dims = 1)]
end
Flux.@functor Fourier

x = rand(Float32, 5, 10000)

F0 = Fourier0(100)
F = Fourier(100)

using BenchmarkTools
@btime F0(x);
@btime F(x);
