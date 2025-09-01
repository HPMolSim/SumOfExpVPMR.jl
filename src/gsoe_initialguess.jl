function gsoe_initialguess(N::Int) 
       Nbig = BigFloat(N)

       # quadrature points θ = π*(1:2:N-1)'/N
       θ = big(π) .* (BigFloat.(1:2:N-1) ./ Nbig)

       # constant b
       b = BigFloat("0.1194")

       # define a BigFloat imaginary unit
       i_big = Complex{BigFloat}(0, 1)

       # contour points z = N*(0.1309 – b*θ.^2 + 0.2500i*θ)
       z = Nbig .* (BigFloat("0.1309") .- b .* θ.^2 .+ BigFloat("0.2500") .* i_big .* θ)

       # derivatives zp = N*(-2b*θ + 0.2500i)
       zp = Nbig .* (-2b .* θ .+ BigFloat("0.2500") .* i_big)

       # weights cs = (1i/N)*exp(z).*zp
       cs = (i_big ./ Nbig) .* exp.(z) .* zp
       
       w = -sqrt(big(π)) .* cs ./ sqrt.(z)
       s = 2*sqrt.(z)

       s = vcat(s, conj.(s))
       w = vcat(w, conj.(w))

       return s, w
end


function gsoe_initialguess1(N::Int) 
       Nbig = BigFloat(N)

       # quadrature points θ = π*(1:2:N-1)'/N
       θ = big(π) .* (BigFloat.(1:2:N-1) ./ Nbig)

       # constant b
       b = BigFloat("0.1194")

       # define a BigFloat imaginary unit
       i_big = Complex{BigFloat}(0, 1)

       # contour points z = N*(0.1309 – b*θ.^2 + 0.2500i*θ)
       z = Nbig .* (BigFloat("0.1309") .- b .* θ.^2 .+ BigFloat("0.2500") .* i_big .* θ)

       # derivatives zp = N*(-2b*θ + 0.2500i)
       zp = Nbig .* (-2b .* θ .+ BigFloat("0.2500") .* i_big)

       # weights cs = (1i/N)*exp(z).*zp
       cs = (i_big ./ Nbig) .* exp.(z) .* zp
       
       w = -sqrt(big(π)) .* cs ./ sqrt.(z)
       s = 2*sqrt.(z)

       s = vcat(s, conj.(s))
       w = vcat(w, conj.(w))

       return s, w
end

function gsoe_initialguess2(N::Int) 
       Nbig = BigFloat(N)

       # quadrature points θ= π*(1:2:N-1)'/N
       θ = big(π) .* (BigFloat.(1:2:N-1) ./ Nbig)

       # define a BigFloat imaginary unit
       i_big = Complex{BigFloat}(0, 1)

       z  = BigFloat("2.246")*Nbig*(BigFloat("1.0") .- sin.(BigFloat("1.1721") .- BigFloat("0.3443").*i_big.*θ));
       zp = BigFloat("0.3443").*i_big.*BigFloat("2.246").*Nbig*cos.(BigFloat("1.1721").-BigFloat("0.3443").*i_big.*θ);


       # weights cs = (1i/N)*exp(z).*zp
       cs = (i_big ./ Nbig) .* exp.(z) .* zp
       
       w = -sqrt(big(π)) .* cs ./ sqrt.(z)
       s = 2*sqrt.(z)

       s = vcat(s, conj.(s))
       w = vcat(w, conj.(w))

       return s, w
end
