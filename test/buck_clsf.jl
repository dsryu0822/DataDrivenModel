using Flux

function sinm(m,x)
    return [sin(k*x) for k in 1:m]
end
function cosm(m,x)
    return [cos(k*x) for k in 1:m]
end
function fourierDense(dims)
    In, Out = dims
    m = (Out ÷ 2In)
    function f(x)
        @assert length(x) == In throw(DimensionMismatch)
        fourierBasis = vcat(x, [sinm.(m,x)..., cosm.(m,x)...]...)
        if Out < length(fourierBasis)
            return fourierBasis[1:Out]
        else
            return [fourierBasis; zeros(Out - length(fourierBasis))]
        end
    end
end


