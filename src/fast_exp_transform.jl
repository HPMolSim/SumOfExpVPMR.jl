
"""
function FGT1d(q_1, q_2, x, soepara, α)

Calculate the summation sum_{i,j} q_1[i] * q_2[j] * f(|x[i] - x[j]| / α) in O(N) steps, via SOE approximation.

"""
function FET1d(q_1::Vector{TQ}, q_2::Vector{TQ}, x::Vector{T}, soepara::SoePara{TC}; α::T = one(T), sort_x::Vector{Int64} = sortperm(x)) where{TQ, T, TC}

    sum_result = zero(ComplexF64)

    N = length(x)
    @assert length(q_1) == length(q_2) == N

    for i in 1:N
        sum_result += q_1[i] * q_2[i]
    end

    for (s, w) in soepara.sw
        A = ComplexF64(q_2[sort_x[N]])
        sum_result += w * q_1[sort_x[N - 1]] * exp(s * (x[sort_x[N - 1]] - x[sort_x[N]]) / α) * A
        for i in N-2:-1:1
            l = sort_x[i]
            m = sort_x[i + 1]
            A = A * exp( - s * (x[sort_x[i + 2]] - x[sort_x[i + 1]]) / α) + q_2[sort_x[i + 1]]
            sum_result += w * q_1[l] * exp(s * (x[l] - x[m]) / α) * A
        end

        B = ComplexF64(q_2[sort_x[1]])
        sum_result += w * q_1[sort_x[2]] * exp(- s * (x[sort_x[1]] - x[sort_x[2]]) / α) * B
        for i in 3:N
            l = sort_x[i]
            m = sort_x[i - 1]
            B = B * exp(s * (x[sort_x[i - 2]] - x[sort_x[i - 1]]) / α) + q_2[sort_x[i - 1]]
            sum_result += w * q_1[l] * exp(s * (x[m] - x[l]) / α) * B
        end
    end

    return sum_result
end