include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5;

@time _trng = factory_soft(DataFrame, 0.1, tspan = [0, 10], ic = [0, .0446272, -0.119564]; dt)
# _trng.gt = 3 .- 2(_trng.u .> 0.05) .- (_trng.u .< -0.05)
_trng.gt = 1 .+ (abs.(_trng.u) .> 0.05)
θ_ = logrange(1e-20, 1e+2, 11)
h_ = round.(Int64, logrange(1, 10000, 21))
# _trng = _trng[10:10:end, :]

nss__ = []
for h = round.(Int64, logrange(1, 10000, 21))
    trng = deepcopy(_trng[h:h:end, :])
    nss_ = []
    @showprogress for θ = θ_
        try
            add_subsystem!(trng, vrbl, cnfg; θ)
            push!(nss_, length(unique(trng.subsystem)))
        catch
            @error "error in $(θ)"
            push!(nss_, 0)
        end
        # ranking = sort(combine(groupby(trng, :subsystem), nrow), :nrow, rev = true)
        # replace!(trng.subsystem, (ranking.subsystem .=> [sort(unique(ranking.subsystem)); 0; 0][1:2])...)
        # acc = count(trng.gt .== trng.subsystem) / nrow(trng)
        # push!(acc_, acc)
    end
    push!(nss__, nss_)
end
# broadcast.(last, nss__)
scatter(dt*logrange(1, 10000, 11), nss__[1], xscale = :log10)

CSV.write("subsystems.csv", DataFrame(stack(nss__), :auto), bom = true)

# plot(θ_, acc_, xscale = :log10, yformatter = y -> "$(round(100y, digits = 2))%")
# plot(θ_, nss_, xscale = :log10, ylabel = "number of subsystems", legend = :none, xticks = [1e-20, 1e+2], xlabel = "θ")
# CSV.write("labeling.csv", DataFrame(θ = θ_, acc = acc_), bom = true)
CSV.write("subsystems.csv", DataFrame(θ = θ_, nss = last.(nss_)), bom = true)

nss_ = []
@showprogress for h = h_
    trng = deepcopy(_trng[h:h:end, :])
    add_subsystem!(trng, vrbl, cnfg; θ)
    push!(nss_, length(setdiff(unique(trng.subsystem), [0])))
end
# broadcast.(last, nss__)
scatter(dt*h_, nss_, xscale = :log10, xlabel = L"h", ylabel = "number of subsystems", legend = :none)
png("number_of_subsystems")