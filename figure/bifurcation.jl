include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
bfcn = CSV.read("bifurcation/soft_bifurcation.csv", DataFrame)
p1 = scatter(bfcn.hrzn, bfcn.vrtc, legend = :none, yticks = [.5, -.5], xticks = [], xlims = [.1, .3], ms = 1, msw = 0, ma = 0.1, color = :black);
bfcn = CSV.read("bifurcation/soft_bifurcation_rcvd.csv", DataFrame)
p2 = scatter(bfcn.hrzn, bfcn.vrtc, legend = :none, yticks = [.5, -.5], xticks = [], xlims = [.1, .3], ms = 1, msw = 0, ma = 0.1, color = 1);

data = CSV.read("lyapunov/soft_lyapunov.csv", DataFrame)
rcvd = CSV.read("lyapunov/soft_lyapunov_rcvd.csv", DataFrame)
p3 = plot(legend = :none, xticks = [.1, .3], xlims = [.1, .3], yticks = [.5, 0, -2.5])
plot!(p3, data.bp, data.λ1, color = :gray)
plot!(p3, data.bp, data.λ2, color = :gray)
plot!(p3, data.bp, data.λ3, color = :gray)
plot!(p3, rcvd.bp, rcvd.λ1, color = 1)
plot!(p3, rcvd.bp, rcvd.λ2, color = 1)
plot!(p3, rcvd.bp, rcvd.λ3, color = 1)

plot(p1, p2, p3, layout = (3, 1), size = (800, 800), dpi = 300);
png("p123")

@info "------------------------------------------------------------------"

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
bfcn = CSV.read("bifurcation/gear_bifurcation.csv", DataFrame)
p1 = scatter(bfcn.hrzn, bfcn.vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, color = :black);
bfcn = CSV.read("bifurcation/gear_bifurcation_rcvd.csv", DataFrame)
p2 = scatter(bfcn.hrzn, bfcn.vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, color = :red);
plot(p1, p2, layout = (2, 1), size = (800, 400))
