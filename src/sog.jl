# 1) a plain, data‐only struct — no inner constructors here
struct SOGPara{T}
    s::Vector{T}
    w::Vector{T}
    function SOGPara(s::Vector{T}, w::Vector{T}) where T
        length(s) == length(w) ||
            throw(ArgumentError("SoePara: length(s) = $(length(s)) ≠ length(w) = $(length(w))"))
        new{T}(s, w)
    end
end

# sum-of-exponentials evaluator
function sog(x::T1, p::SOGPara{T2}; T::DataType = Float64) where {T1<:Real, T2<:AbstractFloat}
    total = zero(T2)
    @inbounds for i in eachindex(p.s)
        total += p.w[i] * exp(-p.s[i] * x.^2)
    end
    return T(real(total))
end

function sog(x::T, s::Vector{T1}, w::Vector{T2}) where{T, T1, T2}
    return sum([w[i] * exp(-s[i] * x.^2) for i in 1:length(s)])
end

function sog_error(f::Function, s::Vector{T2}, w::Vector{T2}, x::Vector{T1} ) where{T1<:Real, T2}
    error = [abs(sog(xᵢ, s, w) - f(xᵢ)) for xᵢ in x]
    return error
end

function sog_max_error(f::Function, s::Vector{T2}, w::Vector{T2}, x::Vector{T1}) where{T1<:Real, T2}
    error = sog_error(f, s, w, x)
    return maximum(error)
end

# error vector for sum‐of‐exponentials approximation
function sog_error(f::Function,
                   p::SOGPara{T},
                   x::Vector{T1}
                  ) where {T1<:Real, T<:AbstractFloat}
    [abs(soe(xᵢ, p) - f(xᵢ)) for xᵢ in x]
end

# maximum error over the same grid
function sog_max_error(f::Function,
                   p::SOGPara{T},
                   x::Vector{T1}
                  ) where {T1<:Real, T<:AbstractFloat}
    maximum(sog_error(f, p, x))
end
