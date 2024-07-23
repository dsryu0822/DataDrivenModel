vrbl = [:dV, :dI], [:V, :I]
cnfg = (; N = 1)
dt = 1e-7
θ1 = 1e+1; θ2 = 1e+0; θ3 = 1e+0; min_rank = 2;

vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3

vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-2)
dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;