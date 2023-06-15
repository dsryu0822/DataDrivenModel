using LinearAlgebra

function col_normalize(M)
    return M ./ norm.(eachcol(M))'
end