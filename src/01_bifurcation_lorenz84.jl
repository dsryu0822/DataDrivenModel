include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-84

''''''''''''''''''''''''''''''''''''''''''''''''''"""

F0, G0 = 3, 0.2
F1, G1 = 0.05, 0.2 # Hopf bifurcation
F1, G1 = 5, 1.0 # chaotic transition
F1, G1 = 9, 0.5 # period doubling
F_ = range(F0, F1, length = 11)
G_ = range(G0, G1, length = 11)
sol0 = factory_lorenz84(DataFrame, [F_[ 1 ], G_[ 1 ]]; ic = [1, 1, 1], saveat = 900:1e-2:1000)
sol1 = factory_lorenz84(DataFrame, [F_[end], G_[end]]; ic = [1, 1, 1], saveat = 900:1e-2:1000)
plot(
    plot(sol0.x, sol0.y, sol0.z, alpha = .5, color = :black),
    plot(sol1.x, sol1.y, sol1.z, alpha = .5, color = :black),
)

using MATLAB
mat"""
plot3($(sol1.x), $(sol1.y), $(sol1.z))
"""

# bfcn = callbfcn()
traj_ = [DataFrame() for _ in eachindex(F_)]
@showprogress @threads for k in eachindex(F_)
    F, G = F_[k], G_[k]
    sol = factory_lorenz84(DataFrame, [F, G]; ic = [1, 1, 1], saveat = 900:1e-2:1000)
    # bfcn[k] = sol.z[arglmax(sol.z)]
    traj_[k] = sol
end
# scatter(dict2bifurcation(bfcn)..., color = :red, msw = 0, ms = 3, ma = 1)

ys = [traj_[k].y for k in eachindex(F_)]
zs = [traj_[k].z for k in eachindex(F_)]
n  = length(F_)
@mput ys zs n

mat"""
fig = figure('Position', [100 100 800 400]);
hold on;
for k = 1:11
    x = k * ones(10000, 1);
    plot3(x, $ys{k}, $zs{k}, 'k');
end
hold off;
ax = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.ZTickLabel = [];

pbaspect([4 1 1]);   % x축을 y, z축보다 4배 길게 (핵심)
grid on;
view(-15, 15);
"""

mat"""
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);

Z =  3*(1-X).^2.*exp(-(X.^2) - (Y+1).^2) ...
   - 10*(X/5 - X.^3 - Y.^5).*exp(-X.^2-Y.^2) ...
   - 1/3*exp(-(X+1).^2 - Y.^2) ...
   + 4*exp(-((X-2).^2 + (Y-2).^2));   % 새로 추가한 봉우리

surf(X, Y, Z);
shading interp;
colormap(parula);
lighting gouraud;
camlight;
axis off;
view(-30, 45);
"""
