include("../core/header.jl")

function lyapunov_exponent(_data::DataFrame, J_, DT::Root{Float64, Int64}, bf_param;
    U = I(ncol(_data)), T = (last(_data.t) - first(_data.t)))

    λ = zeros(size(U, 1))
    for k = 1:nrow(_data)
        s = apply_tree(DT, collect(_data[k, :]))
        J = J_[s](collect(_data[k, :]))
        U, V = gram_schmidt(U)
        λ += V |> eachcol .|> norm .|> log
        U = RK4(J, U, dt)
    end
    return sort(λ / T, rev=true)
end
wsum(a, b, s = 0) = (1 - s) * a + s * b

function dt2df(dtree)
    io = IOBuffer()
    print_tree(io, dtree)
    str_dtree = replace(String(take!(io)), r"\w.\d+/\d+\w." => "",
                        "Feature " => "", "?" => "", " : " => "", "    "=>"\t", " "=>"")
    lines = split(str_dtree, "\n")[1:(end-1)]
    
    virtualdata = DataFrame([Float64[] for _ in last(vrbl)], last(vrbl))
    virtualdata.subsystem = Int64[]
    
    depth = (length.(findall.("\t", lines)) + length.(findall.("<", lines)))
    depth[1] = 0
    for dpth = sort(unique(depth), rev = true)
        bit_dpth = (depth .== dpth)
        leaf = strip.(lines[bit_dpth], Ref(['\t', '├', '└', '─']))
        featarg = parse(Int64, split(leaf[1], "<")[1])
        featval = parse(Float64, split(leaf[1], "<")[2])
        push!(virtualdata, [[ifelse(featarg == k, featval - 128eps(), 0) for k in eachindex(last(vrbl))]; parse(Int64, leaf[2])])
        push!(virtualdata, [[ifelse(featarg == k, featval + 128eps(), 0) for k in eachindex(last(vrbl))]; parse(Int64, leaf[3])])
        depth[bit_dpth] .= -1
        depth[findlast(bit_dpth)] = dpth - 1
    end
    sort!(virtualdata, last(vrbl))
    return virtualdata
end
function bifurcation_(sysname::AbstractString, data, d)
    if sysname == "soft"
        idx_sampled = diff(abs.(data.u) .> (d/2)) .< 0
        return data[Not(1), :v][idx_sampled]    
    elseif sysname == "gear"
        idx_sampled = diff([0; mod.(data.Ω, 2π)]) .< 0
        return data.v[idx_sampled]
    elseif sysname == "hrnm"
        idx_sampled = idx_sampled = abs.(diff(data.dz)) .> 0.1
        return data[Not(1), :x][idx_sampled]
    end
end


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
sysname = "soft" # d ∈ [0.1, 0.3]
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi], λ = 2e-1) # λ = 5e-1 → 1e-2 → 1e-3
dt = 1e-5; tspan = [30, 50]; θ = 1e-6;


_idx = 1:2001
# done = unique(vcat(CSV.read.(filter(x -> occursin("bfcn", x), readdir()), DataFrame)...).idcs)
# tasks = Dict("chaos1" => 1:650, "chaos2" => 651:1300, "chaos3" => 1301:length(_idx), "SickGPU" => setdiff(_idx, done))
tasks = Dict("chaos1" => 1:10:650, "chaos2" => 651:10:1300, "chaos3" => 1301:10:length(_idx), "Sickbook" => 1:100:length(_idx))
schedules = DataFrame(idx = _idx, bp = LinRange(0.1, 0.3, length(_idx)))[tasks[device], :]
for k in eachindex(last(vrbl)) schedules[!, "λ$k"] .= 0.0 end

data = factory_(sysname)(DataFrame, 0.1; tspan, dt)
add_subsystem!(data, vrbl, cnfg; θ)

_f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
_cnfg_ = [getproperty.(Ref(f), f |> propertynames) for f in _f_];
_Dtree = dryad(data, last(vrbl))
_df_Dtree = dt2df(_Dtree);
_df_Dtree.subsystem .= [3, 1, 1, 2]

# jacobian.(Matrix, f__[1])[3]
J_ = [
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) 0 0],
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) -160000 -172.363],
    vv -> [0 0 0; 0 0 1; -π*sinpi(vv[1]) -160000 -172.363],
]
# M0 = Matrix(_df_Dtree[:, 1:(end-1)])
@info "$(now()): Preprocess for $(sysname) done!"

idcs = [Int64[] for _ in _idx]
hrzn = [Float64[] for _ in _idx]
vrtc = [Float64[] for _ in _idx]
@showprogress @threads for dr = eachrow(schedules)
# try
    # pin = (dr.bp - bp1) / (bp2 - bp1)

    M0 = deepcopy(_df_Dtree)
    M0.u .= [-(dr.bp/2) - 128eps(), -(dr.bp/2) + 128eps(), +(dr.bp/2) - 128eps(), +(dr.bp/2) + 128eps()]
    # M0.subsystem = df_Dtree_[2].subsystem
    Dtree = dryad(M0, last(vrbl)) # print_tree(Dtree)

    f_ = deepcopy(_f_)
    f_[2].sparse_matrix[1, 3] = 80000*dr.bp
    f_[3].sparse_matrix[1, 3] = -80000*dr.bp
    f_[2].dense_matrix[1, 3] = 80000*dr.bp
    f_[3].dense_matrix[1, 3] = -80000*dr.bp

    # J_ = []; while true try J_ = jacobian.(Function, f_); break; catch; print("."); end end
    data = DataFrame(solve(f_, [eps(), .05853, .47898], dt, 0:dt:150, Dtree), last(vrbl));
    data = data[(nrow(data) ÷ 5):end, :]

    λ = lyapunov_exponent(data[:, last(vrbl)], J_, Dtree, dr.bp)
    dr[names(schedules)[3:end]] .= λ
    
    vrtc[dr.idx] = bifurcation_(sysname, data, dr.bp)
    hrzn[dr.idx] = repeat([dr.bp], length(vrtc[dr.idx]))
    idcs[dr.idx] = repeat([dr.idx], length(vrtc[dr.idx]))

    try
        CSV.write("...$(device)ing $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
        bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
        CSV.write("...$(device)ing $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)
    catch e
        @error "\n$(now()): $e in $(dr.idx)"
    end
end
CSV.write("!$(device) $(sysname)_lyapunov_rcvd.csv", schedules, bom = true)
bfcn = DataFrame(idcs = vcat(idcs...), hrzn = vcat(hrzn...), vrtc = vcat(vrtc...))
CSV.write("!$(device) $(sysname)_bfcn_rcvd.csv", bfcn, bom = true)

sysname = "gear" # Fe ∈ [-0.2, 0.2]
sysname = "hrnm" # l ∈ [0, 1]
# plot(data.u[1:100:end])