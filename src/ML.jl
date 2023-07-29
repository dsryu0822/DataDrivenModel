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

