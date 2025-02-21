solve(_f_)


hhh_ = factory_soft(DataFrame, 0.1, ic = [0, .000129715, .301469]; tspan = [0, 12], dt)
@time hhh4 = DataFrame(solve(_f_, collect(hhh_[1, 1:3]), (10dt), 0:(10dt):12, _Dtree), last(vrbl));
@time hhh5 = DataFrame(solve(_f_, collect(hhh_[1, 1:3]), (dt), 0:(dt):12, _Dtree), last(vrbl));
@time hhh6 = DataFrame(solve(_f_, collect(hhh_[1, 1:3]), (dt/10), 0:(dt/10):12, _Dtree), last(vrbl));

plot(
    plot(hhh_.t[1:10:end], hhh_.u[1:10:end], title = "(original) h = 10^(-5)", color = :black, legend = :none),
    plot(hhh4.t[1:1:end], hhh4.u[1:1:end], title = "h = 10^(-4)", color = :black, legend = :none),
    plot(hhh5.t[1:10:end], hhh5.u[1:10:end], title = "h = 10^(-5)", color = :black, legend = :none),
    plot(hhh6.t[1:100:end], hhh6.u[1:100:end], title = "h = 10^(-6)", color = :black, legend = :none),
    layout = (:, 1), size = [800, 800]
)
png("hhhhhhhhhhhh.png")
plot(hhh_.t[1:10:end], hhh_.u[1:10:end], title = "(original) h = 10^(-5)", color = :black, legend = :none)
plot!(hhh5.t[1:10:end], hhh5.u[1:10:end], title = "h = 10^(-5)", color = :blue, legend = :none)