struct GaussParameter{T}
    sw::Vector{NTuple{2, T}}
end

function GaussParameter(Step::Int)  
    legendre_sw = legendre(Step)
    return GaussParameter{Float64}([tuple(legendre_sw[1][i], legendre_sw[2][i]) for i in 1:Step])
end

function GaussParameter(T::DataType, Step::Int) 
    legendre_sw = legendre(T, Step)
    return GaussParameter{T}([tuple(legendre_sw[1][i], legendre_sw[2][i]) for i in 1:Step])
end

@inline function Gauss_int(integrand::Function, Gaussian::GaussParameter{TG}; region::NTuple{2, T} = ((-1.0), (1.0))) where {T <: Number, TG <: Number}

    a, b = TG.(region)

    result = zero(TG)
    for i in 1:size(Gaussian.sw, 1)
        result += integrand((b + a) / TG(2) + (b - a) * Gaussian.sw[i][1] / TG(2)) * (b - a) * Gaussian.sw[i][2] / TG(2)
    end
    
    return result
end

@inline function Gauss_int_vector(integrand::Function, Gaussian::GaussParameter{TG}; region::NTuple{2, T} = (-1.0, 1.0)) where {T <: Number, TG <: Number}

    a, b = TG.(region)

    result = zeros(TG, size(Gaussian.sw, 1))
    for i in 1:size(Gaussian.sw, 1)
        result[i] = integrand((b + a) / TG(2) + (b - a) * Gaussian.sw[i][1] / TG(2)) * (b - a) * Gaussian.sw[i][2] / TG(2)
    end
    
    return result
end