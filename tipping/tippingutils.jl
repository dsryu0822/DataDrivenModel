
function syntheticSINDy(Ξ, sysms::Tuple, recipe::AbstractDataFrame)
    Ysyms, Xsyms = sysms
    bit_sparse = all.(map(x -> iszero.(x), eachrow(Ξ)))
    recipeF = recipe[.!bit_sparse, :]
    _Ξ = Ξ[.!bit_sparse, :]
    mse = 0
    aic = -Inf
    r2 = 1
    return STLSQresult("synthetic", recipe, recipeF, Ξ, _Ξ, mse, aic, r2, Ysyms, Xsyms)
end
function arglmax(x)
    bits = circshift(x, 1) .< x .> circshift(x, -1)
    bits[1] = false
    bits[end] = false
    return findall(bits)
end