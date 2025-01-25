include("../core/header.jl")

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
dt = 1e-5 / 10; tspan = [30, 50];

@time trng = factory_soft(DataFrame, 0.1; tspan, dt)
trng.ss = Int64.(abs.(trng.u) .> 0.05)

##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
dt = 1e-3 / 100; tspan = [0, 1000];

@time trng = factory_hrnm(DataFrame, 0.1; tspan, dt)
trng.ss = Int64.(abs.(trng.z) .> 1)

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
dt = 1e-2 / 100; tspan = [0, 1000];

@time trng = factory_gear(DataFrame, -0.2; tspan, dt)
trng.ss = Int64.(abs.(trng.x) .> 1)

# plot(trng.x[1:100:end], color = trng.ss[1:100:end])
sp_ = unique(round.(Int64, logrange(1, 10000, 50)))
recall = []
for sparsity = sp_
    smpl = trng[1:sparsity:end, :]
    smpl.hit = (smpl.ss .!= circshift(smpl.ss, 1)) + (smpl.ss .!= circshift(smpl.ss, -1))
    push!(recall, count(findall(smpl.hit .== 1) .∈ Ref(detect_jump(smpl, vrbl, dos = 1))))
end
num_discontinuous = sum((trng.ss .!= circshift(trng.ss, 1)) + (trng.ss .!= circshift(trng.ss, -1)))
scatter(dt .* sp_, recall ./ num_discontinuous, color = :black, lw = 1, xscale = :log10, shape = :x, alpha = .5, legend = :none, xlabel = "h", ylims = [-.05, 1.05], xticks = exp10.(-10:10))
png("partition")
