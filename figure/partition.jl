include("../core/header.jl")


##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################

vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3; θ = 1e-10

@time trng = factory_hrnm(DataFrame, 0.1, tspan = [0, 100]; dt)
CSV.write("G:/DDM/partition.csv", trng)
findall(abs.(diff(trng.dz)) .> 1e-1)

temp = findall(abs.(diff(trng.dz)) .> 1e-1)
ground_truth = sort([temp .+ 1; temp .- 1; temp])

using Plots
default(label = :none, color = :black, xlims = [0, nrow(trng)])
pargs = (; xticks = [], yticks = [])
pz = plot(trng.z; color = :lightgray, pargs..., yticks = [-1, 1])
scatter!(pz, ground_truth, trng.z[ground_truth], shape = :x, ms = 2)
dz_withNaN = deepcopy(trng.dz); dz_withNaN[ground_truth] .= NaN
pdz = plot(1:nrow(trng),dz_withNaN; msw = 0, ms = .5, color = :lightgray, pargs...);
pdz_gt = scatter(pdz, ground_truth, trng.dz[ground_truth], shape = :x, ms = 2);

#############################
using L1TrendFiltering
@time x = l1tf(trng.dz, 100).x
# plot(abs.(diff(diff(x))))
l1tf_bps = findall(abs.(diff(diff(x))) .> 1e-3) .+ 1
pdz_l1tf = scatter(pdz, l1tf_bps, trng.dz[l1tf_bps], shape = :x, ms = 2, color = 1);
pdz_l1tf
# pdz_l1tf = plot(pdz_l1tf, x, ls = :dash, color = 1, alpha = 0.5)


#############################
# BPs_dim_1 = trunc.(Int64, 1000CSV.read("C:/Users/rmsms/OneDrive/lab/piecewise regression/BPs_dim=1.csv", DataFrame, header = 0)[2:(end-1), 1])
BPs_dim_2 = trunc.(Int64, 1000CSV.read("C:/Users/rmsms/OneDrive/lab/piecewise regression/BPs_dim=2.csv", DataFrame, header = 0)[2:(end-1), 1])

# pdz_dim_1 = scatter(pdz, BPs_dim_1, trng.dz[BPs_dim_1], shape = :x, ms = 2, color = 2)
pdz_dim_2 = scatter(pdz, BPs_dim_2, trng.dz[BPs_dim_2], shape = :x, ms = 2, color = 2)

#############################
@time result = bm(trng, vrbl, cnfg)[2:end-1]
pdz_proposed = scatter(pdz, result, trng.dz[result], shape = :x, ms = 2, color = 3);
normeddf = [0; norm.(eachrow(diff(diff(Matrix(trng[:, first(vrbl)]), dims = 1), dims = 1))); 0]
pddz = scatter(log10.(normeddf), msw = 0, ms = 1, color = 3; pargs...);

#############################
pparition = plot(pdz_gt, pdz_l1tf, pdz_dim_2, pdz_proposed, pddz, layout = (:, 1), size = [600, 600], dpi = 300);
png(pparition, "partition")

intersect(ground_truth, l1tf_bps)
intersect(ground_truth, BPs_dim_2)
intersect(ground_truth, result)

