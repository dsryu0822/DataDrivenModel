include("../core/header.jl")

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
# function J_(t, u, v, d)
#     return [   0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
# end
schedules = CSV.read("schedules/soft.csv", DataFrame)
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi], λ = 2e-1) # λ = 5e-1 → 1e-2 → 1e-3
dt = 1e-5; tspan = [30, 50]; θ = 1e-6;
@showprogress @threads for dr = eachrow(schedules)
    try
        filename1 = "data/soft/$(lpad(dr.idx, 5, '0')).csv"
        if !isfile(filename1)
            state = factory_soft(DataFrame, dr.bp; tspan, dt)[:, last(vrbl)]
            M_state = Matrix(state)
            derivative = (M_state[2:(end-0),:] - M_state[1:(end-1),:]) ./ dt
            data = [state[1:(end-1), :] DataFrame(derivative, first(vrbl))]
            add_subsystem!(data, vrbl, cnfg; θ)
            if 0 ∉ data.subsystem
                CSV.write(filename1, data, bom = true)

                filename2 = replace(filename1, "data/soft" => "output/soft")
                if !isfile(filename2)
                    f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
                    Dtree = dryad(data, last(vrbl))
                    rcvd = solve(f_, [data[1,1:3]...], dt, 0:dt:abs(-(tspan...)), Dtree)
                    CSV.write(filename2, DataFrame(rcvd, last(vrbl)), bom = true)
                end
            else
                continue
            end
        end
    catch
        open("error.csv", "a") do io println(io, "$(now()),$(dr.idx)") end
        @error "\n$(now()): error in $(dr.idx)"
    end
end
@showprogress @threads for dr = eachrow(schedules)
    try
        filename1 = "data/soft/$(lpad(dr.idx, 5, '0')).csv"
        filename2 = replace(filename1, "data/soft" => "output/soft")
        data = CSV.read(filename1, DataFrame)
        f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
        Dtree = dryad(data, last(vrbl))
        rcvd = solve(f_, [data[1,1:3]...], dt, 0:dt:abs(-(tspan...)), Dtree)
        CSV.write(filename2, DataFrame(rcvd, last(vrbl)), bom = true)
    catch
        open("error.csv", "a") do io println(io, "$(now()),$(dr.idx)") end
        @error "\n$(now()): error in $(dr.idx)"
    end
end

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
# function J_(x, v, Ω, θ, Fe)
#     dfdx = ifelse(abs(x) > 1, 1, 0)
#     return [ 0 1 0 0
#              -(1 + k1*cos(Ω))*dfdx -2ζ (-Fe*H*sin(Ω) + k1*sin(Ω)*dfdx) -Fe*H*sin(θ)
#              0 0 0 0
#              0 0 0 0 ]
# end
schedules = CSV.read("schedules/gear.csv", DataFrame)
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-2)
# λ = 1e-2 works for 762 bps/ λ = 1e-4 works for 138 bps
dt = 1e-2; tspan = [0, 1000]; θ = 1e-3; dos = 1

idx_tgt = Not(parse.(Int64, first.(readdir("output/gear"), 5)))
schedules = CSV.read("schedules/gear.csv", DataFrame)[idx_tgt, :]
# dr = eachrow(schedules)[1]
@showprogress @threads for dr = eachrow(schedules)
    try
        filename1 = "data/gear/$(lpad(dr.idx, 5, '0')).csv"
        if !isfile(filename1)
            # _data = factory_gear(DataFrame, dr.bp; tspan, dt)
            state = factory_gear(DataFrame, dr.bp; tspan, dt)[:, last(vrbl)]
            M_state = Matrix(state)
            derivative = (M_state[2:(end-0),:] - M_state[1:(end-1),:]) ./ dt
            data = [state[1:(end-1), :] DataFrame(derivative, first(vrbl))][20000:end, :]
            add_subsystem!(data, vrbl, cnfg; θ, dos)
            if 0 ∉ data.subsystem
                filename2 = replace(filename1, "data/gear" => "output/gear")
                if !isfile(filename2)
                    f_ = [SINDy(df[rand(1:end, 20), :], vrbl...; cnfg...) for df in groupby(data, :subsystem)]
                    Dtree = dryad(data, last(vrbl))
                    rcvd = solve(f_, [data[1,1:4]...], dt, 0:dt:abs(-(tspan...)), Dtree)
                    CSV.write(filename1, data, bom = true)
                    CSV.write(filename2, DataFrame(rcvd, last(vrbl)), bom = true)
                end
            else
                print("$(dr.idx), ")
                continue
            end
        end
    catch
        open("error.csv", "a") do io println(io, "$(now()),$(dr.idx)") end
        @error "\n$(now()): error in $(dr.idx)"
    end
end
# f_ = [SINDy(df[rand(1:end, 20), :], vrbl...; cnfg...) for df in groupby(data, :subsystem)]
# f = SINDy(data[data.subsystem .== 1, :][rand(1:end, 20),:], vrbl...; cnfg...)
# print(f)
# print.(f_)
# plot(data.x, color = data.subsystem)
# _data = factory_gear(DataFrame, dr.bp; tspan, dt)[20000:end, :]
# data[findall(diff(_data.x .> 1) .!= 0), :]
# data.subsystem |> unique
# plot(data.x, color = data.subsystem, xlims = [0, 10000])
# scatter!(jumpt, data.x[jumpt])
# data[data.subsystem .== 2, :]
# subsystem |> unique


CSV.write("241226/diff.csv", DataFrame(diff = sum(abs2, diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 2)[1:100000]))
CSV.write("241226/dat.csv", data)
plot(sum(abs2, diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 2)[1:100000])