using Plots, LaTeXStrings
@info "Packages Plots, LaTeXStrings loaded"
default(msw=0, color=:black)
mm = Plots.mm

function look(data)
    vn = names(data)
    _data = data[1:100:end, :]
    p2 = plot(_data.t, _data[:,2]   , xlabel =  L"t", ylabel = vn[2], legend = :none)
    p3 = plot(_data.t, _data[:,3]   , xlabel =  L"t", ylabel = vn[3], legend = :none)
    p4 = plot(_data[:,2], _data[:,3], xlabel = vn[2], ylabel = vn[3], legend = :none)
    return plot(p3, p4, plot(), p2, size = (800, 800))
end