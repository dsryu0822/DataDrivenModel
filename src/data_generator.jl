include("../core/header.jl")

function lyapunov_exponent(_data, J_, bf_param;
    U = Matrix(qr(J_(collect(_data[1, :])..., bf_param)).Q))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        J = J_(collect(_data[k, :])..., bf_param)
        U, V = gram_schmidt(U)
        if 2k > nrow(_data)
            λ += V |> eachcol .|> norm .|> log
        end
        U = RK4(J, U, dt)
    end
    return sort(2λ / (last(data.t) - first(data.t)))
end

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
schedules = CSV.read("bifurcation/soft_schedules.csv", DataFrame)
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2)
dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
function J_(t, u, v, d)
    return [   0                                0                                 0
               0                                0                                 1
     -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
end
@showprogress @threads for dr = eachrow(schedules)
    filename = "bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)
    # data = factory_soft(DataFrame, dr.idx, dr.d, tspan = [0, 50]); data = data[30(nrow(data) ÷ 50):end , :]
    add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
    CSV.write(filename, data)

    # idx_sampled = diff(abs.(data.u) .> (dr.d/2)) .> 0
    # sampledv = data[Not(1), :v][idx_sampled]
    # append!(hrzn, fill(dr.d, length(sampledv)))
    # append!(vrtc, sampledv)
end

##########################################################################
#                                                                        #
#                           DC-DC buck converter                         #
#                                                                        #
##########################################################################
schedules = CSV.read("bifurcation/buck_schedules.csv", DataFrame)
vrbl = [:dV, :dI], [:V, :I]
cnfg = (; N = 1)
dt = 1e-7; θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;
function J_(V, I, Vr, E)
    _R = 22 # 22
    _L = 20_m # 20m
    _C = 47_μ # 22μ
    return [ -1/(_R*_C) (1/_C)
                -(1/_L)      0 ]
end

@showprogress @threads for dr = eachrow(schedules)
    filename = "bifurcation/buck/$(lpad(dr.idx, 5, '0')).csv"
    # data = CSV.read(filename, DataFrame)
    data = factory_buck(DataFrame, dr.idx, dr.E, tspan = [0, .30]); # data = data[29(nrow(data) ÷ 30):end , :]
    add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank); # 30 sec
    CSV.write(filename, data)

    # idx_sampled = diff(data.Vr) .< 0
    # sampledV = data[Not(1), :V][idx_sampled]
    # append!(hrzn, fill(dr.E, length(sampledV)))
    # append!(vrtc, sampledV)
end

##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
# schedules = CSV.read("bifurcation/hrnm_schedules.csv", DataFrame)
schedules = CSV.read("lyapunov/hrnm_schedules_cache.csv", DataFrame)
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3; θ1 = 1e-2; θ2 = 1e-27; θ3 = 1e-1; min_rank = 32;
function J_(t, x, y, z, _f)
     _a,  _b,  _c,  _d,  _k,  _ω,  _α,  _β = 
    1.0, 3.0, 1.0, 5.0, 0.9, 1.0, 0.1, 0.8
    return [                0                             0   0    0
             -_ω*_f*sin(_ω*t) (-3*_a*(x^2) + 2*_b*x + _k*z)   1 _k*x
                            0                       -2*_d*x  -1    0
                            0                            _β   0  -_α]
end

schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; schedules[!, :λ4] .= .0;


@showprogress @threads for dr = eachrow(schedules)[.!isfile.(["lyapunov/hrnm_traj/$(lpad(dr.idx, 5, '0')).csv" for dr = eachrow(schedules)])]
    try
        # filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
        filename = "lyapunov/hrnm_traj/$(lpad(dr.idx, 5, '0')).csv"
        # data = CSV.read(filename, DataFrame)
        # data = factory_hrnm(DataFrame, dr.idx, dr.f, tspan = [0, 1500]); data = data[1000(nrow(data) ÷ 1500):end , :]
        data = factory_hrnm(DataFrame, dr.idx, dr.f, ic = [dr.t, dr.x, dr.y, dr.z], tspan = [0, 500]; dt = 1e-4)
        add_subsystem!(data, vrbl, cnfg; θ1, θ2, θ3, min_rank)
        # CSV.write(filename, data)
        CSV.write(replace(filename, "bifurcation/hrnm" => "lyapunov/hrnm_traj"), data)

        # idx_sampled = abs.(diff(data.dz)) .> 0.1
        # sampledx = data[Not(1), :x][idx_sampled]
        # append!(hrzn, fill(dr.f, length(sampledx)))
        # append!(vrtc, sampledx)
        λ = lyapunov_exponent(data[:, last(vrbl)], J_, dr.f)
        dr[[:λ1, :λ2, :λ3, :λ4]] .= λ
    catch
        @error "Error in $(dr.idx)"
    end
end
CSV.write("lyapunov/hrnm.csv", schedules, bom = true)
# plot(data.z[1:100:end], color = data.subsystem[1:100:end])




# schedules[!, :t] .= .0; schedules[!, :x] .= .0; schedules[!, :y] .= .0; schedules[!, :z] .= .0;
# for dr = eachrow(schedules)
#     filename = "bifurcation/hrnm/$(lpad(dr.idx, 5, '0')).csv"
#     data = CSV.read(filename, DataFrame)
#     dr[(end-3):end] .= collect(data[end, 1:4])
# end
# CSV.write("lyapunov/hrnm_schedules_cache.csv", schedules, bom = true)