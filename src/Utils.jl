using FiniteDifferences

function fdiff(x::T; stencil = 5, dt = 1, method = central_fdm) where T <: AbstractVector
    fdm = method(stencil, 1)
    dx = reduce(+, circshift.(Ref(x), -fdm.grid) .* fdm.coefs) / dt
    for boundary in fdm.grid
        if boundary < 0
            popfirst!(dx)
        elseif boundary > 0
            pop!(dx)
        end
    end
    return dx
end
function fdiff(X::T; stencil = 5, dt = 1, method = central_fdm) where T <: AbstractMatrix
    fdm = method(stencil, 1)
    DX = vcat(transpose.(fdiff.(eachrow(X); stencil = stencil, dt = dt, method = method))...)
    for boundary in fdm.grid
        if boundary < 0
            X = X[:, 1:(end-1)]
        elseif boundary > 0
            X = X[:, 2:end]
        end
    end
    return DX, X
end

function Euler(f::Function,v::Array{Float64,1}, h=10^(-2))
    return v + h*f(v)
end
function RK4(f::Function,v::Array{Float64,1}, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4)
end
