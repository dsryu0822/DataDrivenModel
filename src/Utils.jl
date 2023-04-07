using FiniteDifferences

function fdiff(x::T; order = 5, dt = 1) where T <: AbstractVector
    fdm = backward_fdm(order, 1)
    dx = reduce(+, circshift.(Ref(x), -fdm.grid) .* fdm.coefs) / dt
    for boundary in fdm.grid
        if boundary > 0
            pop!(dx)
        elseif boundary < 0
            popfirst!(dx)
        end
    end
    return dx
end

function fdiff(X::T; order = 5, dt = 1) where T <: AbstractMatrix
    # return vcat(transpose.(fdiff.(Array{Float64}.(eachrow(X))))...)
    return vcat(transpose.(fdiff.(eachrow(X); order = order, dt = dt))...)
end


backward_fdm(2, 1)
# using Plots
rand(10)'
transpose(X)
fdiff(X)
fdiff(X[1, 1:10])
dt = 0.1
x = 0:dt:2π
y = sin.(x)
dy1 = fdm.(sin, x)[4:(end-3)]
dy2 = fdiff(y, order = 7, dt = dt)

plot((dy1 - dy2), label = "error")
# reduce(+ , circshift.(Ref(X), 0, fdm.grid) .* fdm.coefs)
