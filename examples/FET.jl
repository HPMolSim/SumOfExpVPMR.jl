using Plots, SumOfExpVPMR

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

    # warm up
    direct_sum(f, rand(2), rand(2), rand(2))
    FET1d(rand(2), rand(2), rand(2), soepara)


    N_array = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600]
    for N in N_array
        q_1 = rand(N)
        q_2 = rand(N)
        x = rand(N)

        sort_x = sortperm(x)

        push!(time_direct, @elapsed direct_sum(f, q_1, q_2, x))
        push!(time_FET, @elapsed FET1d(q_1, q_2, x, soepara, sort_x = sort_x))
    end
end

begin
    plot(xlabel = "log(N)", ylabel = "log(time)", dpi = 500, title = "Compare direct sum and FET")
    plot!(log10.(N_array), log10.(time_direct), label = "direct sum", marker = :circle)
    plot!(log10.(N_array), log10.(time_FET), label = "FET", marker = :square)
    savefig("FET.png")
end