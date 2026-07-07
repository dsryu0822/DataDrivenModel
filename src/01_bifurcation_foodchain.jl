include.("../core/" .* readdir("core")[[1,2,3,4,6]])

multi_bifurcation = plot()
for jld2 in filter(x ->  endswith(x, ".jld2"), readdir())
    @info "Loading bifurcation data for food chain from file..."
    bfcn = JLD2.load(jld2)["bfcn"]
    scatter!(multi_bifurcation, dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, xticks = [.9, .961, 1.0])
    scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, xticks = [.9, .961, 1.0])
    png("G:/food/$(jld2)_scatter.png")
end
png(multi_bifurcation, "G:/food/multi_bifurcation_scatter.png")

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
K_ = 0.88:1e-4:1.00
c0_ = 12e-2:1e-2:46e-2
@showprogress @threads for c0 in c0_
    if !isfile("bifurcation_foodchain_$(rpad(c0, 4, '0')).jld2")
        @info "Calculating bifurcation data for food chain with c0 = $(c0)..."
        bfcn = Dict{Float64, Vector{Float64}}()
        for k in eachindex(K_)
            sol = factory_foodchain(DataFrame, K_[k], ic = [.85, c0, .8], saveat = 0:1e-1:10000)
            P_ = sol.P[sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[K_[k]] = P_[arglmin(P_)]
            end
        end
        JLD2.@save "bifurcation_foodchain_$(rpad(c0, 4, '0')).jld2" bfcn
    else
        @info "Loading bifurcation data for food chain from file..."
        JLD2.@load "bifurcation_foodchain_$(rpad(c0, 4, '0')).jld2" bfcn
    end
end

sol = factory_foodchain(DataFrame, 0.961, ic = [.85, 0.10, .8], saveat = 0:1e-1:10000)[(end-1000):end, [:R, :C, :P]]
plot(sol.R, sol.C, sol.P, lims = [0, 1.1])
CSV.write("sol10.csv", sol)

arglmin(sol2.P)
collect(sol2[102, :])

@showprogress for K in 0.9600:1e-6:0.9620
    sol1 = factory_foodchain(DataFrame, K, ic = [0.7440555124967513, 0.2666843335892154, 0.7601257636315785], saveat = 0:1e-1:10000)[(end-1000):end, [:R, :C, :P]]
    sol2 = factory_foodchain(DataFrame, K, ic = [0.7470781897589035, 0.26629728666884384, 0.7381698314917428], saveat = 0:1e-1:10000)[(end-1000):end, [:R, :C, :P]]
    plot(title = "K = $(K)", size = [400, 400], dpi = 300)
    plot!(sol1.R, sol1.C, sol1.P, color = :blue, lw = 2, alpha = 0.5)
    plot!(sol2.R, sol2.C, sol2.P, color = :red, lw = 2, alpha = 0.5)
    png("C:/Users/rmsms/Downloads/새 폴더/foodchain_K_$(rpad(K, 8, '0')).png")
end
gif(anim, "foodchain_bifurcation.gif", fps = 10)
CSV.write("sol10.csv", sol)