shrinkage(x::Real, α::Real) = sign(x) * max(abs(x)-α, zero(x))

function soft_thresholding(x, λ)
    return sign(x) * max(abs(x) - λ, 0.0)
end

function ista(A, b, λ = 1e-1, α = 1e-12, max_iter=1000)
    L₂ = norm.(eachcol(A))
    A = A ./ L₂'

    x = A \ b
    for i in 1:max_iter
        grad = A' * (A * x - b)
        x = soft_thresholding.(x - α * grad, α * λ)
    end
    return sparse(x ./ L₂)
end
data = factory_lorenz(DataFrame, 28, tspan = [0, 100])

A = Θ(data[1:100000:end, [:x, :y, :z]]; N = 2)

rank(A)
ista(A, data[1:100000:end, :dx], 1e-1, 1e-12, 1000)

A = Θ(datasets[24][:, last(vrbl)]; B_[6]...)
ista(A, datasets[24].dt, 1e-5, 1e-1, 1000)

@info ""
# TODO fast iterative shrinkage and hard-thresholding algorithm
# for weighted l1-norm minimization
function fista(A::AbstractMatrix, b::AbstractVector, λ::Real,
    x::AbstractVector = spzeros(size(A, 2)); maxiter::Int = 1024, stepsize::Real = 1e-2)
    w = fill(λ, size(x))
    fista(A, b, w, x, maxiter = maxiter, stepsize = stepsize)
end

function l1(x::AbstractVector, w::AbstractVector)
    length(x) == length(w) || throw(DimensionMismatch("length(x) ≠ length(w)"))
    n = zero(eltype(x))
    @simd for i in eachindex(x)
        @inbounds n += w[i] * abs(x[i])
    end
    return n
end

# TODO: stepsize selection
function ista(A, b, λ::Real, x = spzeros(size(A, 2)); maxiter::Int = 1024, stepsize::Real = 1e-2)
    ista(A, b, fill(λ, size(x)), x, maxiter = maxiter, stepsize = stepsize)
end

function ista(A::AbstractMatrix, b::AbstractVector, w::AbstractVector,
                                        x::AbstractVector = spzeros(size(A, 2));
                                        maxiter::Int = 1024, stepsize::Real = 1e-2)
    x = sparse(x)
    r(x) = b-A*x # residual
    f(x) = sum(abs2, r(x)) + l1(x, w)
    g(x) = A'r(x) # negative gradient
    α = stepsize
    fx = f(x)
    for i in 1:maxiter
        ∇ = g(x)
        @. x = shrinkage(x + 2α*∇, w*α)
        dropzeros!(x) # optimize sparse representation
    end
    return x
end


function fista(A::AbstractMatrix, b::AbstractVector, w::AbstractVector,
    x::AbstractVector = spzeros(size(A, 2));
    maxiter::Int = 1024, stepsize::Real = 1e-2)
    x = sparse(x)
    y = copy(x)
    t = one(eltype(x))
    for i in 1:maxiter
        ∇ = A' * (b - A * y)
        x_new = shrinkage.(y .+ 2 * stepsize * ∇, w .* stepsize)
        t_new = (1 + sqrt(1 + 4 * t^2)) / 2
        y = x_new + ((t - 1) / t_new) * (x_new - x)
        x = x_new
        t = t_new
        dropzeros!(x)
    end
    return x
end

cnfg = (; N = 2)
temp = factory_lorenz(DataFrame, 28)
sampled = shuffle(1:nrow(temp))[1:5]
A = Θ(temp[sampled, [:x, :y, :z]]; cnfg...)

fista(A, temp.dy[sampled], 100; stepsize = 1e-8)
