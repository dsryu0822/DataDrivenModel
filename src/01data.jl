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
# idx_tgt = Not(parse.(Int64, first.(readdir("output/soft"), 5)))
# schedules = CSV.read("schedules/soft.csv", DataFrame)[idx_tgt, :]
schedules = CSV.read("schedules/soft.csv", DataFrame)
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2) # λ = 5e-1 → 1e-2 → 1e-3
dt = 1e-5; tspan = [30, 50]; θ = 1e-6;

# dr = eachrow(schedules)[1159]
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
                CSV.write(filename1, data)
            else
                continue
            end
        else
            data = CSV.read(filename1, DataFrame)
        end

        filename2 = replace(filename1, "data/soft" => "output/soft")
        if !isfile(filename2)
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
            Dtree = dryad(data, last(vrbl))
            rcvd = solve(f_, [data[1,1:3]...], dt, 0:dt:20, Dtree)
            CSV.write(filename2, data)
        end
    catch
        open("error.csv", "a") do io
            println(io, "$(now()),$(dr.idx)")
        end
        @error "$(now()): error in $(dr.idx)"
    end
end
# unique(data.subsystem)
# values(check_)

# for filename in ["output/soft/$(lpad(idx, 5, '0')).csv" for idx in [keys(check_)...][0 .∈ values(check_)]]
#     try
#         rm(filename)
#     catch
#         @error "$(now()),$(filename),remove"
#     end
# end

# setdiff(readdir("G:/DDM/data/soft"), readdir("G:/DDM/output/soft"))
plot(data.u[1:100:end], color = data.subsystem[1:100:end])
plot(plot(data.u[200000:201100], yticks = [-dr.bp, dr.bp]/2), plot(data.v[200000:201100]), layout = (2,1))
