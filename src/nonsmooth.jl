function argjump(type::Type, Y::AbstractVector)
    Δ²Y = diff(diff(Y))
    Δ²Y = Δ²Y ./ maximum(Δ²Y)
    bit_jump = Δ²Y .> 0.1
    println(sum(bit_jump), " jumps detected!")
    if type == Bool
        return bit_jump
    elseif type == Int64
        return findall(bit_jump)
    else
        @error "type must be Bool or Int64"
    end
end
argjump(Y) = argjump(Int64, Y) # If type is not specified, return indices