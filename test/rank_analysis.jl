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
    temp = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
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
data = factory_hrnm(DataFrame, 0.1)
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]

cnfg = forwardselect(cook(last(vrbl), poly = 0:4, trig = 0:2, format = sin), data, vrbl)
labeling!(data, vrbl, cnfg, θ = 1e-5)
f_ = [SINDy(subdata, vrbl, cnfg, λ = 1e-3) for subdata in groupby(data, :label)]
f_ .|> print
