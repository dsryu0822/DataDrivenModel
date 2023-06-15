using Flux
using DataFrames

struct Affine
    W
    b
end
Affine(in::Integer, out::Integer) = Affine(randn(out, in), randn(out))
(m::Affine)(x) = m.W * x .+ m.b
Flux.@functor Affine
a = Affine(10, 5)
a(rand(10))
Flux.params(a)
