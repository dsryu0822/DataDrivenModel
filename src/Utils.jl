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