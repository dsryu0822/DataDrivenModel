
totex = f0

replace(
    replace(
        define(totex, sigdigits = 5),
        reverse(totex.recipe.tex[2:end] .=> totex.recipe.term[2:end])...,
    ), "⋅" => "", "²" => "^{2}", "³" => "^{3}", "⁴" => "^{4}", "⁵" => "^{5}",
) |> print