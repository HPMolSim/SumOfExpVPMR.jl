# 1) a plain, data‐only struct — no inner constructors here
struct SoePara{T}
    s::Vector{T}
    w::Vector{T}
    function SoePara(s::Vector{T}, w::Vector{T}) where T
        length(s) == length(w) ||
            throw(ArgumentError("SoePara: length(s) = $(length(s)) ≠ length(w) = $(length(w))"))
        new{T}(s, w)
    end
end

# sum-of-exponentials evaluator
function soe(x::T1, p::SoePara{T2}; T::DataType = Float64) where {T1<:Real, T2}
    xabs = abs(x)
    total = zero(T2)
    @inbounds for i in eachindex(p.s)
        total += p.w[i] * exp(-p.s[i] * xabs)
    end
    return T(real(total))
end

function soe(x::T, s::Vector{T1}, w::Vector{T2}) where{T, T1, T2}
    return sum([w[i] * exp(-s[i] * x) for i in 1:length(s)])
end

function soe_error(f::Function, s::Vector{T2}, w::Vector{T2}; x::Vector{T1} = big.([0.0:0.01:10.0...])) where{T1<:Real, T2}
    error = [abs(soe(xᵢ, s, w) - f(xᵢ)) for xᵢ in x]
    return error
end

function max_error(f::Function, s::Vector{T2}, w::Vector{T2}; x::Vector{T1} = big.([0.0:0.01:10.0...])) where{T1<:Real, T2}
    error = soe_error(f, s, w, x = x)
    return maximum(error)
end

# error vector for sum‐of‐exponentials approximation
function soe_error(f::Function,
                   p::SoePara{T};
                   x::Vector{T1} = big.([0.0:0.01:10.0...])
                  ) where {T1<:Real, T}
    [abs(soe(xᵢ, p) - f(xᵢ)) for xᵢ in x]
end

# maximum error over the same grid
function max_error(f::Function,
                   p::SoePara{T};
                   x::Vector{T1} = big.([0.0:0.01:10.0...])
                  ) where {T1, T}
    return maximum(soe_error(f, p; x = x))
end
