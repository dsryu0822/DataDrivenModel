using LinearAlgebra, Clustering, CSV, DataFrames

function col_normalize(M)
    return M ./ norm.(eachcol(M))'
end

struct Fourier
    m::Int64
    an::AbstractArray{Float32, 1}
    bn::AbstractArray{Float32, 1}
    a0::AbstractArray{Float32, 1}
    L::AbstractArray{Float32, 1}
end
function Base.show(io::IO, l::Fourier)
    print(io, "Fourier(", l.m, ")")
end
Fourier(m::Integer) = Fourier(m, zeros(Float32, m), zeros(Float32, m), zeros(Float32, 1), zeros(Float32, 1))
function (layer::Fourier)(x)
    # a0_L = exp.(layer.a0_L)
    nxL⁻¹ = (1:layer.m) .* x[[end], :] .* exp.(layer.L)

    # return [x[1:end-1, :]; layer.a0 .+ sum(
    return [x; layer.a0 .+ sum(
        (cospi.(nxL⁻¹) .* layer.an) + 
        (sinpi.(nxL⁻¹) .* layer.bn)
        , dims = 1)]
end
Flux.@functor Fourier

struct Modulo
    T::AbstractArray{Float32, 1}
end
function Base.show(io::IO, l::Modulo)
    print(io, "Modulo(", l.T, ")")
end
Modulo(T::AbstractFloat) = Modulo(Float32[T])

function (layer::Modulo)(x)
    return [x[1:(end-1),:]; mod.(x[[end],:], layer.T)]
end
Flux.@functor Modulo



isgpu(M::AbstractArray) = M isa CuArray
isgpu(f::Chain) = any([typeof(first(f)).types...] .<: CuArray)
gcpu(M) = ifelse(isgpu(M), gpu, cpu)

function Accuracy(f, data, label)
    x, _ = data.data
    n = size(x, 2)
    m = length(label)
    @assert n == m "data($n) and label($m) size mismatch"

    x = x |> gcpu(f)
    predicted = argmax.(eachcol(f(x) |> cpu))
    return sum(label .== predicted) / n
end

function init_ANN(data)
    p = size(data.data[1], 1)
    ANN = Chain( # SubSystemSelector
        # Flux.Scale(p, bias = false),
        # Fourier(100),
        # Dense(p+1=> 50, relu),
        # Modulo(400*(10^(-6))),
        Dense(p => 100, relu),
        Dense(100 => 100, relu),
        Dense(100 => 100, relu),
        Dense(100 => nsubsys),
        softmax
        ) |> gcpu(data.data[1])
    @info "ANN is initiating on $(ifelse(isgpu(ANN), "GPU", "CPU"))"

    Loss(x,y) = Flux.crossentropy(ANN(x), y)
    return ANN
end
function save_ANN(ANN, data; lastepch = 0, dir = "C:/Temp/")
    @info "ANN is training on $(ifelse(isgpu(ANN), "GPU", "CPU"))"
    optimizer = ADAM()
    ps = Flux.params(ANN)

    x, y = data.data
    Loss(x, y) = Flux.crossentropy(ANN(x), y)

    loss_ = [Loss(x, y)]
    acry_ = [Accuracy(ANN, data, label)]
    epch_ = [0]

try
    @info "Training started!"
    for epch in 1:lastepch
        if (epch < 10) || (epch % 1000 == 0)
            @time Flux.train!(Loss, ps, data, optimizer)
        else
            Flux.train!(Loss, ps, data, optimizer)
        end
        loss = Loss(x, y)
        if loss < minimum(loss_)
            acry = Accuracy(ANN, data, label)
            if acry > maximum(acry_)
                @info join(["epoch ", lpad(epch, 5)
                , ": loss = ", rpad(trunc(loss, digits = 6), 8)
                , ", acry = ", rpad(trunc(100acry, digits = 4), 8)
                , ","])
    
                print("\033[F!\r\n")
                push!(loss_, loss)
                push!(acry_, acry)
                push!(epch_, epch)
                _ANN = ANN |> cpu
                jldsave(joinpath(dir, "SSS_$(lpad(epch, 6, '0')).jld2"); _ANN, loss, acry, epch)
                jldsave(joinpath(dir, "SSS.jld2"); _ANN, loss, acry, epch)
            end
        end
    end
    cp("C:/Temp/SSS.jld2", "data/SSS.jld2")
catch err
    if err isa InterruptException
        @info """
        Stopped by you. the last one is backed up:
        $(dir)/_ANN.jld2
        """
    else
        @error err
    end
    jldsave(joinpath(dir, "_ANN.jld2"); _ANN = ANN, loss = loss_[end], acry = acry_[end], epch = epch_[end])
end
    return ANN
end
function load_ANN(dir)
    file = jldopen(dir)
    ANN = file["_ANN"]
    loss = file["loss"]
    acry = file["acry"]
    epch = file["epch"]
    close(file)
    @info """
    ANN is loaded on $(ifelse(isgpu(ANN), "GPU", "CPU"))
    loss = $loss
    acry = $acry
    epch = $epch
    """
    return ANN, loss, acry, epch
end
