include("../core/header.jl")

##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
schedules = CSV.read("schedules/soft.csv", DataFrame)
idcs, hrzn, vrtc = Int64[], Float64[], Float64[]
@showprogress for dr = eachrow(schedules)
    filename = "output/soft/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)

    idx_sampled = diff(abs.(data.u) .> (dr.bp/2)) .> 0
    sampledv = data[Not(1), :v][idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledv)))
    append!(hrzn, fill(dr.bp, length(sampledv)))
    append!(vrtc, sampledv)
end
bfcn = DataFrame(; idcs, hrzn, vrtc)
CSV.write("output/bfcn_soft.csv", bfcn, bom = true)
bfcn = CSV.read("output/bfcn_soft.csv", DataFrame)

scatter(bfcn.hrzn, bfcn.vrtc, msw = 0, ms = .5, color = :black, ylims = [-1, 1]);
png("temp")

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
schedules = CSV.read("schedules/gear.csv", DataFrame)
idcs, hrzn, vrtc = Int64[], Float64[], Float64[]
@showprogress for dr = eachrow(schedules)
    try
    filename = "output/gear/$(lpad(dr.idx, 5, '0')).csv"
    data = CSV.read(filename, DataFrame)

    # data = data[data.Ω .> 5000, :]
    idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
    sampledx = data.v[idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledx)))
    append!(hrzn, fill(dr.bp, length(sampledx)))
    append!(vrtc, sampledx)
    catch
        continue
    end
end
bfcn = DataFrame(; idcs, hrzn, vrtc)
CSV.write("output/bfcn_gear.csv", bfcn, bom = true)
bfcn = CSV.read("output/bfcn_gear.csv", DataFrame)

scatter(bfcn.hrzn, bfcn.vrtc, msw = 0, ms = .5, color = :black, ylims = [-2, 2]);
png("temp")