
bfcn = CSV.read("bifurcation/gear_bifurcation_rcvd.csv", DataFrame)
scatter(bfcn.hrzn, bfcn.vrtc, ms = 1, legend = :none, msw = 0, ma = 0.1, color = :black);
png("bifurcation/gear_bifurcation_rcvd.png")