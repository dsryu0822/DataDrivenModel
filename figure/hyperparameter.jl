include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5; tspan = [0, 10]; θ = 1e-16;
# dt = 1e-5 is good to bifurcation diagram but not for lyapunov spectrum, 1e-6 is needed
# θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
else
    @time data = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
end
# add_subsystem!(data, vrbl, cnfg; θ)
# Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)

# measure1 = []
# measure2 = []
# λ_ = logrange(1e-20, 1e+2, 100)
# @showprogress for λ = λ_
#     cnfg = (; f_ = [cospi], λ = λ)
#     f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)] # print.(f_)
#     push!(measure1, f_[2].MSE)
#     push!(measure2, length(f_[2].sparse_matrix.nzval))
# end

# plot(
#     scatter(λ_, measure1, xscale = :log10, legend = :none, xlabel = L"\lambda", ylabel = "MSE"),
#     scatter(λ_, measure2, xscale = :log10, legend = :none, xlabel = L"\lambda", ylabel = "nzval"),
#     layout = (2, 1), size = [600, 600], dpi = 300
# )
# png("lambda vs 2 things")



# sys_recovered = []
# sys_parameter = []
# performance = []
# test_ = []
# @showprogress for (n, m) = Base.product(0:3, 0:3)
#     # cnfg = (; N = n, M = m, f_ = [cospi], λ = 1e-1)
#     cnfg = (; N = n, M = m, λ = 1e-1)

#     f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)] # print.(f_)
#     test = DataFrame(solve(f_, collect(data[1, last(vrbl)]), dt, first(tspan):dt:last(tspan), Dtree)[1:(end-1),:], last(vrbl))
#     push!(test_, test)

#     push!(performance, sum(abs2, data.u - test.u))
#     push!(sys_recovered, f_)
#     push!(sys_parameter, "(N = $n, M = $m)")
# end
# sys_recovered = reshape(sys_recovered, 4, 4)
# sys_parameter = reshape(sys_parameter, 4, 4)
# performance = reshape(performance, 4, 4)
# # performance = broadcast(x -> sum(getproperty.(x, :MSE)), sys_recovered)
# # nzvals = broadcast(x -> sum(length.(getproperty.(getproperty.(x, :sparse_matrix), :nzval))), sys_recovered)


# p1 = plot(data.u[1:10:end], color = :black, legend = :none)
# for (k, test) in enumerate(test_)
#     if performance[k] < 30
#         plot!(p1, test.u[1:10:end])
#     end
# end
# p1

# sys_recovered[2, 1]
# sys_recovered[end, end]
# log10.(performance)

# print(sys_recovered[1, 3][3])

λ_ = logrange(1e-10, 1e-1, 4)
N_ = 0:3
M_ = 0:3

schedule = DataFrame(ID = Int64[], λ = Float64[], N = Int64[], M = Int64[], ss = Int64[], error = Float64[], runtime = [])
for (n, m, λ) in Base.product(N_, M_, λ_) push!(schedule, (0, λ, n, m, 0, 0.0, 0)) end
schedule.ID = 1:size(schedule, 1)

@showprogress @threads for dr = eachrow(schedule)
    try
        tic = now()
        λ = dr.λ; n = dr.N; m = dr.M
        cnfg = (; N = n, M = m, λ = λ)
        
        trng = deepcopy(data)
        add_subsystem!(trng, vrbl, cnfg; θ)
        if 0 ∈ trng.subsystem
            dr.ss = -1
        else
            dr.ss = length(unique(trng.subsystem))
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :subsystem)] # print.(f_)
            Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
            test = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, first(tspan):dt:last(tspan), Dtree)[1:(end-1),:], last(vrbl))
            dr.error = sum(abs2, trng.u - test.u)
            CSV.write("hyperparameter/$ID n=$(n)_m=$(m).csv", test, bom = true)
        end
        print(schedule)
        CSV.write("65 result.csv", schedule, bom = true)
        toc = now()
        dr.runtime = (now() - tic).value / 60000
    catch
        dr.runtime = -1
    end
end
