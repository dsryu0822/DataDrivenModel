include("../core/header.jl")

doublerange(n) = Base.product(1:n, 1:n)
# has_neighbor(arr) = !isempty(intersect(arr, [arr .+ 1; arr .- 1]))

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
    test = CSV.read("hyperparameter/0 test.csv", DataFrame)
else
    temp = factory_["soft"](DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    data = temp[1:(nrow(temp) ÷ 2), :]
    test = temp[(nrow(temp) ÷ 2):end, :]
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
    CSV.write("hyperparameter/0 test.csv", test, bom = true)
end
vrbl = [:dt, :du, :dv], [:t, :u, :v]

cnfg = forwardselect(cook(last(vrbl), poly = 0:3, trig = 0:2, format = cospi), data, vrbl)
labeling!(data, vrbl, cnfg, θ = 1e+0)
f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-3) for subdata in groupby(data, :label)]
f_ .|> print

##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
data = factory_["hrnm"](DataFrame, 0.1)
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]

cnfg = forwardselect(cook(last(vrbl), poly = 0:4, trig = 0:2, format = sin), data, vrbl)
labeling!(data, vrbl, cnfg, θ = 1e-5)
f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-3) for subdata in groupby(data, :label)]
f_ .|> print

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
    data = CSV.read("buck_r/$(lpad(dr.idx, 5, '0')).csv", DataFrame)

    bit_12 = data.V .< 11.958
    bit_poincare = (bit_12 .!= circshift(bit_12, 1)) .&& (data.I .< 0.55)
    push!(vtcl, data.I[bit_poincare])
    push!(hrzn, fill(dr.bp, sum(bit_poincare)))
end
scatter([hrzn...;], [vtcl...;], msw = 0, ms = 1, color = :black, legend = :none)


data = factory_["buck"](DataFrame, 40)
vrbl = [:dV, :dI], [:V, :I]

cnfg = forwardselect(cook(last(vrbl), poly = 0:2, trig = 0:2), data, vrbl)
labeling!(data, vrbl, cnfg, θ = 1e-3)
f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-9) for subdata in groupby(data, :label)]
f_ .|> print
# atree = dryad(data, vrbl)
@showprogress @threads for dr = eachrow(schedules)[1:10:end]
    E = dr.bp
    property_ = deepcopy([getproperty(f_[2], property) for property in propertynames(f_[2])])
    property_[3][1, 2] = 50*E
    property_[4][1, 2] = 50*E
    _f_ = [f_[1], STLSQresult(property_...)]

    rcvd = DataFrame(solve(_f_, collect(data[1, last(vrbl)]), 0:1e-7:0.5, atree), last(vrbl))
    # rcvd = DataFrame(solve(_f_, [11.5, .9], 0:1e-7:0.5, atree), last(vrbl))
    rcvd = rcvd[4(nrow(rcvd) ÷ 5):end, :]
    # CSV.write("buck_r/$(lpad(dr.idx, 5, '0')).csv", rcvd)
    # rcvd = DataFrame(solve(_f_, [11.5, 0.9], 0:1e-7:0.5, atree), last(vrbl))
    # rcvd = rcvd[4(nrow(rcvd) ÷ 5):end, :]
    CSV.write("buck_re/$(lpad(dr.idx, 5, '0')).csv", rcvd)
end
@showprogress for dr = eachrow(schedules)[1:10:end]
    E = dr.bp
    rcvd = CSV.read("buck_re/$(lpad(dr.idx, 5, '0')).csv", DataFrame)
    plot(rcvd.V, rcvd.I, size = [400, 400], legend = :none, color = 1, xlims = [11.5, 14], ylims = [.2, .9], title = E)
    # plot!([11.958, 11.958], [minimum(rcvd.I), 0.55], color = 2)
    png("buck_re/$(lpad(dr.idx, 5, '0')).png")
end

idx = 2
property_ = [getproperty(f_[2], property) for property in propertynames(f_[2])]
property_[3][1, 2] = 50*eachrow(schedules)[idx].bp
property_[4][1, 2] = 50*eachrow(schedules)[idx].bp
_f_ = [f_[1], STLSQresult(property_...)]
_f_ .|> print

rcvd = DataFrame(solve(_f_, collect(data[1, last(vrbl)]), 0:1e-7:0.5, atree), last(vrbl))
rcvd = rcvd[4(nrow(rcvd) ÷ 5):end, :]
plot(rcvd.V, rcvd.I, size = [400, 400], legend = :none, color = 1)
rcvd = CSV.read("buck_r/$(lpad(idx, 5, '0')).csv", DataFrame)
plot!(rcvd.V, rcvd.I, size = [400, 400], legend = :none, color = 2)

DecisionTree.confusion_matrix(prune_tree(atree, 0.9))
DecisionTree.confusion_matrix(data.label, apply_tree(prune_tree(atree, 0.6), Matrix(data[:, [:V, :I]])))


data = factory_["buck"](DataFrame, 40, tspan = [0, 1])
vrbl = [:dV, :dI], [:V, :I]
cnfg = cook(last(vrbl), poly = 0:1)
# labeling!(data, vrbl, cnfg, θ = 1e-3)
data[!, :label] .= (data.V .< data.Vr) .+ 1
f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-9) for subdata in groupby(data, :label)]
f_ .|> print
_data = data[1:100000, :]
_data[(_data.label .!= circshift(_data.label, -1)) .|| (_data.label .!= circshift(_data.label, 1)), :]
@time atree = dryad(_data, vrbl)

# rcvd = DataFrame(solve(f_, [11.5, .9], 0:1e-7:0.5, atree), last(vrbl))
@time rcvd = DataFrame(solve(f_, collect(data[1, last(vrbl)]), 0:1e-7:0.5, atree), last(vrbl))
plot(rcvd.V[4000000:10:end], rcvd.I[4000000:10:end], size = [400, 400], legend = :none, color = 1, label = :none)
scatter!(rcvd.V[[4000000]], rcvd.I[[4000000]], legend = :bottomright, label = "t = 0.4")
scatter!(rcvd.V[[end]], rcvd.I[[end]], label = "t = 0.5")

accuracy = nfoldCV_tree(data.label, Matrix(data[:, last(vrbl)]), 5)
3.2*10*log(10)

f_ .|> print
_f_ .|> print
f_[2].matrix[1,2]
_f_[2].matrix[1,2]

build_forest(data.label, Matrix(data[:, last(vrbl)]), rng = 1)

plot(rcvd.V[0000001:10:2000000], rcvd.I[0000001:10:2000000], size = [400, 400])
plot(rcvd.V[2000001:10:4000000], rcvd.I[2000001:10:4000000], size = [400, 400])
plot(rcvd.V[4000001:10:5000000], rcvd.I[4000001:10:5000000], size = [400, 400])


data = factory_["buck"](DataFrame, 40, ic = [11.5, .7], tspan = [0, 0.005])
qck1 = plot(data.V[1:10:end], color = :black, label = "V", xlims = [0, 5000])
plot!(qck1, data.Vr[1:10:end], color = :red, label = "Vr")
qck2 = plot(data.dI[1:10:end], color = :black, label = "dI", xlims = [0, 5000])
plot(qck1, qck2, layout = (2, 1), plot_title = "original system")

# quick = DataFrame(solve(f_, collect(data[1, last(vrbl)]), 0:1e-7:0.005, atree), last(vrbl))
quick = DataFrame(solve(f_, [11.5, .7], 0:1e-7:0.005, atree), last(vrbl))
qck1 = plot(quick.V[1:10:end], color = :black, label = "V", xlims = [0, 5000])
plot!(qck1, data.Vr[1:10:end], color = :red, label = "Vr")
qck2 = plot(diff(quick.I)[1:10:end], color = :black, label = "dI", xlims = [0, 5000])
plot(qck1, qck2, layout = (2, 1), plot_title = "decision tree with 10^5")

foo = dryad(data, vrbl; method = :forest)
typeof(foo) <: Root
typeof(foo) <: Ensemble
solve(f_, collect(data[1, last(vrbl)]), 0:1e-7:0.005, foo)