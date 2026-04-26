using Pkg
Pkg.activate("ad_env")
Pkg.instantiate()
using DataFrames
using CairoMakie
using Statistics

kalibracja = DataFrame([[0.259, 0.23, 0.214, 0.28, 0.468, 0.424, 0.219], [0.438, 1.136, 0.455, 0.373, 0.946, 0.561, 0.61], [0.39, 0.228, 0.439, 0.444, 0.208, 0.278, 0.466], [0.226, 0.382, 0.218, 0.266, 0.32, 0.334, 0.385]], [:g1_s1, :g1_s2, :g2_s1, :g2_s2])

test = DataFrame([[0.335, 0.291, 0.379], [1.059, missing, missing], [0.588, 0.362, 0.393], [0.262, 0.252, 0.278]], [:g1_s1, :g1_s2, :g2_s1, :g2_s2])

for c in names(kalibracja)
    test[!, c] .-= kalibracja[1, c]
    kalibracja[!, c] .-= kalibracja[1, c]
end
kalibracja[!, :g1_średnia_wiersze] = mean(Matrix(kalibracja[!, [:g1_s1, :g1_s2]]), dims=2)[:, 1]
kalibracja[!, :g2_średnia_wiersze] = mean(Matrix(kalibracja[!, [:g2_s1, :g2_s2]]), dims=2)[:, 1]
kalibracja[!, :g1g2_średnia_wiersze] = mean(Matrix(kalibracja[!, [:g2_s1, :g2_s2, :g1_s1, :g1_s2]]), dims=2)[:, 1]
kalibracja[!, :próbka] = [0, 1, 2, 4, 6, 8, 10]

function regresja_liniowa(x, y)
    x_śr = mean(x)
    y_śr = mean(y)
    a = sum((x .- x_śr) .* (y .- y_śr)) / sum((x .- x_śr) .^ 2)
    b = y_śr - a * x_śr
    r2 = sum((y_śr .- (a .* x .+ b)) .^ 2) / sum((y_śr .- y) .^ 2)
    return [a, b, r2]
end

function wykres(ax, y, x, test_x)
    mask = ismissing.(x)
    mask = .!mask
    x = x[mask]
    y = y[mask]
    a, b, r2 = regresja_liniowa(x, y)
    linspace = [minimum(skipmissing(vcat(x, test_x))), maximum(skipmissing(vcat(x, test_x)))]
    test_y = test_x .* a .+ b
    vlines!(ax, test_x, linewidth=1.2, linestyle=:dash, color=:red, alpha=0.7)
    lines!(ax, linspace, a .* linspace .+ b, linestyle=(:dash, :dense), label="r²=" * string(round(r2, digits=3)), linewidth=4)
    scatter!(ax, x, y, markersize=13, color=:lightblue)
    scatter!(ax, test_x, test_y, markersize=17, marker=:x, color=:red)
    axislegend(ax, position=:rb, framevisible=false)
end

bkg_col = :gray95

fig = Figure(size=(1200, 600))
ax = Axis(fig[1, 1], title="grupa1 sesja1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], ylabel="próbka")
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g1_s1], test[:, :g1_s1])
ax = Axis(fig[1, 2], title="grupa1 sesja2", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka])
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g1_s2], test[:, :g1_s2])
ax = Axis(fig[2, 1], title="grupa2 sesja1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], ylabel="próbka", xlabel="absorbancja względna")
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g2_s1], test[:, :g2_s1])
ax = Axis(fig[2, 2], title="grupa2 sesja2", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], xlabel="absorbancja względna")
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g2_s2], test[:, :g2_s2])
ax = Axis(fig[1, 3], title="grupa1 średnia", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka])
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g1_średnia_wiersze], vcat(test[:, :g1_s1], test[:, :g1_s2]))
ax = Axis(fig[2, 3], title="grupa2 średnia", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], xlabel="absorbancja względna")
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g2_średnia_wiersze], vcat(test[:, :g2_s1], test[:, :g2_s2]))
#display(fig)
save("plot1.png", fig)

fig = Figure()
ax = Axis(fig[1, 1], title="Średnie wartości pomiarów dla krzywej kalibracyjnej we wszystkich seriach.\nUwzględniono również średnie z testowych dla grupy2 oraz grupy1 serii1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], xlabel="absorbancja względna", ylabel="próbka")
wykres(ax, kalibracja[:, :próbka], kalibracja[:, :g1g2_średnia_wiersze], [mean(test[:, :g2_s1]), mean(test[:, :g2_s2]), mean(test[:, :g1_s1])])
#display(fig)
save("plot2.png", fig)

#=
optymalizacja_r2, odrzucanie najmniej dopasowanego punktu
=#

function optymalizacja_r2(x, y)
    Int64(maximum(reverse!.([pushfirst!(regresja_liniowa(x[1:end.!=i], y[1:end.!=i]), i) for i in 1:length(x)]))[end])
end

bez_punktów = DataFrame([[i != optymalizacja_r2(col, kalibracja[:, :próbka]) ? col[i] : missing for i in 1:length(col)] for col in eachcol(kalibracja[:, 1:4])], names(kalibracja[:, 1:4]))
fig = Figure(size=(800, 800))
ax = Axis(fig[1, 1], title="grupa1 seria1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], ylabel="próbka")
wykres(ax, kalibracja[:, :próbka], bez_punktów[:, :g1_s1], test[:, :g1_s1])
ax = Axis(fig[1, 2], title="grupa1 seria2", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka])
wykres(ax, kalibracja[:, :próbka], bez_punktów[:, :g1_s2], test[:, :g1_s2])
ax = Axis(fig[2, 1], title="grupa2 seria1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], xlabel="absorbancja względna", ylabel="próbka")
wykres(ax, kalibracja[:, :próbka], bez_punktów[:, :g2_s1], test[:, :g2_s1])
ax = Axis(fig[2, 2], title="grupa2 seria2", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], xlabel="absorbancja względna")
wykres(ax, kalibracja[:, :próbka], bez_punktów[:, :g2_s2], test[:, :g2_s2])
#display(fig)
save("plot3.png", fig)


bez_punktów[!, :g1g2_średnia_wiersze] = [mean(skipmissing(v)) for v in eachrow(bez_punktów)]
fig = Figure()
ax = Axis(fig[1, 1], title="Średnie wartości pomiarów dla krzywej kalibracyjnej we wszystkich seriach.\nUwzględniono również średnie z testowych dla grupy2 oraz grupy1 serii1", backgroundcolor=bkg_col, yticks=kalibracja[:, :próbka], ylabel="próbka")
wykres(ax, kalibracja[:, :próbka], bez_punktów[:, :g1g2_średnia_wiersze], [mean(test[:, :g2_s1]), mean(test[:, :g2_s2]), mean(test[:, :g1_s1])])
save("plot4.png", fig)
