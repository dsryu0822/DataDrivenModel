using ProgressBars
packages = [:DataFrames, :CSV, :LinearAlgebra, :Clustering,
            :CUDA, :Flux, :JLD2,
            :Plots, :LaTeXStrings]
@time for package in ProgressBar(packages)
    @eval using $(package)
end
println(join(packages, ", "), " loaded!")

include("../src/DDM.jl")
include("../src/ML.jl")
include("../src/nonsmooth.jl")
include("../src/ODEdata.jl")

default(size = (600,600), color = :black, legend = :none)

## Data Load
DATA = CSV.read("G:/buck/buck_000006.csv", DataFrame)
Y = DATA[end-100000:end,[:dV, :dI]] |> Matrix # .|> Float32
X = DATA[end-100000:end,[ :V,  :I]] |> Matrix # .|> Float32
# # XY = [X Y]

# ## Clustering
@showtime dbs = dbscan(col_normalize(Y)', 0.001);
# nsubsys = 2
nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
# # plot(X[:,1], X[:,2], color = dbs.assignments, alpha = 0.5)
# # scatter(Y[:,1], Y[:,2], color = dbs.assignments, alpha = 0.5, msw = 0)
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ_ = [poly_basis(X[s,:], 2)     for s in bit_]
Y_ = [Y[vcat(findall(s)...),:] for s in bit_]
Ξ_ = [STLSQ(Θ_[s], Y_[s], 0.01) for s in 1:nsubsys]

function foo(v)
    v = deepcopy(v)
    s = Int64(pop!(v))
    Θv = poly_basis(v, 2)'
    return vec([(Θv*Ξ_[s]) 0])
end