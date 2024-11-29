include("../core/header.jl")

# ##########################################################################
# #                                                                        #
# #                            Soft impact model                           #
# #                                                                        #
# ##########################################################################
# # function J_(t, u, v, d)
# #     return [   0                                0                                 0
# #                0                                0                                 1
# #      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0) ]
# # end
# schedules = CSV.read("schedules/soft.csv", DataFrame)
# vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 2e-1) # λ = 5e-1 → 1e-2 → 1e-3
# dt = 1e-5; tspan = [30, 50]; θ = 1e-6;

# @showprogress @threads for dr = eachrow(schedules)
#     try
#         filename1 = "data/soft/$(lpad(dr.idx, 5, '0')).csv"
#         if !isfile(filename1)
#             state = factory_soft(DataFrame, dr.bp; tspan, dt)[:, last(vrbl)]
#             M_state = Matrix(state)
#             derivative = (M_state[2:(end-0),:] - M_state[1:(end-1),:]) ./ dt
#             data = [state[1:(end-1), :] DataFrame(derivative, first(vrbl))]
#             add_subsystem!(data, vrbl, cnfg; θ)
#             if 0 ∉ data.subsystem
#                 CSV.write(filename1, data, bom = true)
#             else
#                 continue
#             end
#         else
#             data = CSV.read(filename1, DataFrame)
#         end

#         filename2 = replace(filename1, "data/soft" => "output/soft")
#         if !isfile(filename2)
#             f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
#             Dtree = dryad(data, last(vrbl))
#             rcvd = solve(f_, [data[1,1:3]...], dt, 0:dt:abs(-(tspan...)), Dtree)
#             CSV.write(filename2, data, bom = true)
#         end
#     catch
#         open("error.csv", "a") do io
#             println(io, "$(now()),$(dr.idx)")
#         end
#         @error "$(now()): error in $(dr.idx)"
#     end
# end

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
dt = 1e-2; tspan = [0, 1000]; θ = 1e-4; dos = 1

# idx_tgt = Not(parse.(Int64, first.(readdir("output/gear"), 5)))
# schedules = CSV.read("schedules/gear.csv", DataFrame)[idx_tgt, :]
# dr = eachrow(schedules)[1]
@showprogress @threads for dr = eachrow(schedules)
    try
        filename1 = "data/gear/$(lpad(dr.idx, 5, '0')).csv"
        if !isfile(filename1)
            state = factory_gear(DataFrame, dr.bp; tspan, dt)[:, last(vrbl)]
            M_state = Matrix(state)
            derivative = (M_state[2:(end-0),:] - M_state[1:(end-1),:]) ./ dt
            data = [state[1:(end-1), :] DataFrame(derivative, first(vrbl))]
            add_subsystem!(data, vrbl, cnfg; θ, dos)
            if 0 ∉ data.subsystem
                CSV.write(filename1, data, bom = true)
            else
                continue
            end
        else
            data = CSV.read(filename1, DataFrame)
        end

        filename2 = replace(filename1, "data/gear" => "output/gear")
        if !isfile(filename2)
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
            Dtree = dryad(data, last(vrbl))
            rcvd = solve(f_, [data[1,1:3]...], dt, 0:dt:abs(-(tspan...)), Dtree)
            CSV.write(filename2, data, bom = true)
        end
    catch
        open("error.csv", "a") do io
            println(io, "$(now()),$(dr.idx)")
        end
        @error "$(now()): error in $(dr.idx)"
    end
end

data.subsystem |> unique
plot(data.x, color = data.subsystem)

sample = data[data.subsystem .== 1, :][rand(1:end, 20),:]
SINDy(sample, vrbl...; cnfg...) |> print

gt = (-1 .< data.x .< 1) + 2(data.x .≥ 1) + 3(data.x .≤ -1)
scatter(gt, data.subsystem, legend = :none)

jumpt = detect_jump(data, vrbl, cnfg; dos)
set_divider([1, 4, 5, 6, 10, 20, 25])

sum(abs2, [data[24707, 5:8]...] .- candy([data[24707, 1:4]...])) < 10θ