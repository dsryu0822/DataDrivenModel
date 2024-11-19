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

scatter(bfcn.hrzn, bfcn.vrtc, msw = 0, ms = .5);
png("temp")