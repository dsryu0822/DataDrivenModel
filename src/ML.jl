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

# using Flux
# struct Fourier0
#     an
#     bn
#     L
# end
# function Base.show(io::IO, l::Fourier0)
#     print(io, "Fourier(", length(l.an), ")")
# end

# Fourier0(m::Integer) = Fourier0(randn(Float32, m), randn(Float32, m), Float32[0.0])
# function (layer::Fourier0)(x)
#     m = length(layer.an)
#     nxL⁻¹ = (1:m) * x[end, :]' .* exp(layer.L[1])
#     return [x[1:end-1, :]; sum(
#         [cospi.(nxL⁻¹) .* layer.an
#        ; sinpi.(nxL⁻¹) .* layer.bn], dims = 1)]
# end
# Flux.@functor Fourier0
# ------------------------------

struct Fourier
    m::Int64
    an::AbstractArray{Float32, 1}
    bn::AbstractArray{Float32, 1}
    a0::AbstractArray{Float32, 1}
    L::AbstractArray{Float32, 1}
end
function Base.show(io::IO, l::Fourier)
    print(io, "Fourier(", l.m, ")")
end
Fourier(m::Integer) = Fourier(m, zeros(Float32, m), zeros(Float32, m), zeros(Float32, 1), zeros(Float32, 1))

# function (layer::Fourier)(x)
#     nxL⁻¹ = (1:layer.m) .* x[[end], :] * exp.(layer.a0_L)[end]

#     return [x[1:end-1, :]; identity.(layer.a0_L)[1] .+ sum(
#         (cospi.(nxL⁻¹) .* layer.an) + 
#         (sinpi.(nxL⁻¹) .* layer.bn)
#         , dims = 1)]
# end

function (layer::Fourier)(x)
    # a0_L = exp.(layer.a0_L)
    nxL⁻¹ = (1:layer.m) .* x[[end], :] .* exp.(layer.L)

    # return [x[1:end-1, :]; layer.a0 .+ sum(
    return [x; layer.a0 .+ sum(
        (cospi.(nxL⁻¹) .* layer.an) + 
        (sinpi.(nxL⁻¹) .* layer.bn)
        , dims = 1)]
end
Flux.@functor Fourier

# vcat(sum(CuArray(rand(Float32, 4, size(x, 2))), dims = 1), x)