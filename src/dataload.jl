@time using CSV, DataFrames, Dates

function fillmissing!(data)
    for col in eachcol(data)
        while true
            bit_missing = ismissing.(col)
            if sum(bit_missing) == 0 break end
            col[bit_missing] .= ((circshift(col, -1) + circshift(col, 1)) / 2)[bit_missing]
            bit_missing = ismissing.(col)
            col[bit_missing] .= circshift(col, 1)[bit_missing]
        end
    end
end

function mydata()
    data = DataFrame()
    DATA_ = []
    for fn = readdir("data")
        data
        itemname = Symbol(first(split(fn, " ")))
        if (itemname == :알루미늄) || (itemname == :오렌지) continue end
        push!(DATA_, CSV.read("data/" * fn, DataFrame))
        DATA_[end].날짜 = Date.(DATA_[end].날짜, dateformat"y- m- d")
        select!(DATA_[end], ["날짜", "종가"])
        rename!(DATA_[end], [:t, itemname])
        if eltype(DATA_[end][:, 2]) <: AbstractString
            DATA_[end][:, 2] = replace.(DATA_[end][:, 2], "," => "")
            DATA_[end][!, 2] = parse.(Float64, DATA_[end][:, 2])
        end
        if isempty(data)
            data = deepcopy(DATA_[end])
        else
            leftjoin!(data, DATA_[end], on = :t, makeunique = true)
        end
    end
    data = data[data.t .< Date(2022, 1, 1), :]
    sort!(data, :t)
    fillmissing!(data)
    @assert data == dropmissing(data) # 데이터 무결성 검사
    dropmissing!(data)
    return data
end

data = mydata()
data = data[data.WTI유 .> 0, :]

for col in eachcol(data)[Not(1)]
    col .= (col ./ maximum(col)) #* 2π
end

trng = data[data.t .< Date(2021, 1, 1), :]
test = data[data.t .≥ Date(2020, 12, 31), :]