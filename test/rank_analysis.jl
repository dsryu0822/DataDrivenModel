include("../core/header.jl")

# doublerange(n) = Base.product(1:n, 1:n)
# # has_neighbor(arr) = !isempty(intersect(arr, [arr .+ 1; arr .- 1]))

# ##########################################################################
# #                                                                        #
# #                            Soft impact model                           #
# #                                                                        #
# ##########################################################################
# if isfile("hyperparameter/0 data.csv")
#     data = CSV.read("hyperparameter/0 data.csv", DataFrame)
#     test = CSV.read("hyperparameter/0 test.csv", DataFrame)
# else
#     temp = factory_["soft"](DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
#     data = temp[1:(nrow(temp) ÷ 2), :]
#     test = temp[(nrow(temp) ÷ 2):end, :]
#     CSV.write("hyperparameter/0 data.csv", data, bom = true)
#     CSV.write("hyperparameter/0 test.csv", test, bom = true)
# end
# vrbl = [:dt, :du, :dv], [:t, :u, :v]

# cnfg = forwardselect(cook(last(vrbl), poly = 0:3, trig = 0:2, format = cospi), data, vrbl)
# labeling!(data, vrbl, cnfg, θ = 1e+0)
# f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-3) for subdata in groupby(data, :label)]
# f_ .|> print

# ##########################################################################
# #                                                                        #
# #                           Hindmarsh-Rose model                         #
# #                                                                        #
# ##########################################################################
# data = factory_["hrnm"](DataFrame, 0.1)
# vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]

# cnfg = forwardselect(cook(last(vrbl), poly = 0:4, trig = 0:2, format = sin), data, vrbl)
# labeling!(data, vrbl, cnfg, θ = 1e-5)
# f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-3) for subdata in groupby(data, :label)]
# f_ .|> print

##########################################################################
#                                                                        #
#                           DC-DC Buck Converter                         #
#                                                                        #
##########################################################################
_idx = 1:201
schedules = DataFrame(idx = _idx, bp = LinRange(20, 40, length(_idx)))
# @showprogress @threads for dr in eachrow(schedules)
#     E = dr.bp
#     data = factory_["buck"](DataFrame, E, tspan = [0.0, 0.5])
#     data = data[4(nrow(data) ÷ 5):end, :]
#     CSV.write("buck/$(lpad(dr.idx, 5, '0')).csv", data)
# end

vtcl = []
hrzn = []
@showprogress for dr = eachrow(schedules)
    data = CSV.read("buck_re 1e+6 with t/$(lpad(dr.idx, 5, '0')).csv", DataFrame)

    # bit_12 = data.V .< 11.958
    # bit_poincare = (bit_12 .!= circshift(bit_12, 1)) .&& (data.I .< 0.55)
    bit_poincare = [false; diff(diff(data.I) / 1e-7) .> 10; false]
    # bit_poincare = [(diff(data.dI) .> 10); false]
    push!(vtcl, data.V[bit_poincare])
    push!(hrzn, fill(dr.bp, sum(bit_poincare)))
end
scatter([hrzn...;][1:100:end], [vtcl...;][1:100:end], msw = 0, alpha = .5, ms = 1, color = :black, legend = :none)
scatter([hrzn...;], [vtcl...;], msw = 0, alpha = .5, ms = 1, color = :black, legend = :none)

data = factory_["buck"](DataFrame, 40)
# vrbl = [:dV, :dI], [:V, :I]
vrbl = [:Vr, :dV, :dI], [:t, :V, :I]

# cnfg = forwardselect(cook(last(vrbl), poly = 0:2, trig = 0:2), data, vrbl)
# labeling!(data, vrbl, cnfg, θ = 1e-3)
cnfg = cook(last(vrbl), poly = 0:1)
data[!, :label] = (data.V .< data.Vr) .+ 1

f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-9) for subdata in groupby(data, :label)]
f_ .|> print
# atree = dryad(data, vrbl)

data = factory_["buck"](DataFrame, 40, tspan = [0.0, .1])
data[!, :label] = (data.V .< data.Vr) .+ 1
@time atree = dryad(data, last(vrbl); method = :forest)
@showprogress @threads for dr = eachrow(schedules)[1:1:end]
    E = dr.bp
    property_ = deepcopy([getproperty(f_[2], property) for property in propertynames(f_[2])])
    property_[3][1, end] = 50*E
    property_[4][1, end] = 50*E
    _f_ = [f_[1], STLSQresult(property_...)]

    rcvd = DataFrame(solve(_f_, collect(data[1, last(vrbl)]), 0:1e-7:0.5, atree), last(vrbl))
    rcvd = rcvd[4(nrow(rcvd) ÷ 5):end, :]
    CSV.write("buck_re 1e+6 with t/$(lpad(dr.idx, 5, '0')).csv", rcvd)
end
@showprogress for dr = eachrow(schedules)[1:1:end]
    E = dr.bp
    rcvd = CSV.read("buck_re 1e+6 with t/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    plot(rcvd.V, rcvd.I, size = [400, 400], legend = :none, color = 1, xlims = [11.5, 14], ylims = [.2, .9], title = E)
    png("buck_re 1e+6 with t/$(lpad(dr.idx, 5, '0')).png")
end

# k = 121
# dr = eachrow(schedules)[k]
# foo = CSV.read("buck/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
# bar = CSV.read("buck_re/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
# plot(
#     plot(foo.V[1:10:end], foo.I[1:10:end], title = "original system", color = 2),
#     plot(bar.V[1:10:end], bar.I[1:10:end], title = "recovered system", color = 2),
#     legend = :none, plot_title = "E = $(dr.bp)"
# )
plot(rcvd.V[1:100:end], rcvd.I[1:100:end])
plot(rcvd.V[1:100:250000], rcvd.I[1:100:250000], xlims = [11.5, 13.1], ylims = [.3, .8])
plot(rcvd.V[250001:100:400000], rcvd.I[250001:100:400000], xlims = [11.5, 13.1], ylims = [.3, .8])