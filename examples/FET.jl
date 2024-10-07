using SumOfExpVPMR
using CairoMakie, LaTeXStrings
using BenchmarkTools

function direct_sum(f::Function, q_1::Vector{T}, q_2::Vector{T}, x::Vector{T}) where{T}
    sum_result = zero(T)
    N = length(x)
    for i in 1:N
        for j in 1:N
            sum_result += q_1[i] * q_2[j] * f(abs(x[i] - x[j]))
        end
    end
    return sum_result
end

begin
    f = x -> exp(-x^2)
    soepara, σ = VPMR_cal(f, 6.0, 60, 200, 16)
end

begin
    time_direct = Float64[]
    time_FET = Float64[]

    N_array = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600]
    for N in N_array
        q_1 = rand(N)
        q_2 = rand(N)
        x = rand(N)

        sort_x = sortperm(x)

        push!(time_direct, @belapsed direct_sum($f, $q_1, $q_2, $x))
        push!(time_FET, @belapsed FET1d($q_1, $q_2, $x, $soepara, sort_x = $sort_x))
        @show N, time_direct[end], time_FET[end]
    end
end

begin
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = L"N", ylabel = "time (s)", yscale = log10, xscale = log10)
    scatter!(ax, N_array, time_direct, label = "Direct sum", color = :blue, markersize = 15, marker = :diamond)
    scatter!(ax, N_array, time_FET, label = "FGT", color = :red, markersize = 15, marker = :utriangle)
    x_array = [50:100:25600 * 2...]
    lines!(ax, x_array, x_array.^2 .* (time_direct[end] / N_array[end]^2) , color = :blue, linestyle = :dash)
    lines!(ax, x_array, x_array .* (time_FET[end] / N_array[end]), color = :red, linestyle = :dash)
    axislegend(ax, position = :lt)
    xlims!(ax, 100 * 0.75, 25600 / 0.75)
    fig
    save("FET.png", fig, px_per_unit = 2)
end