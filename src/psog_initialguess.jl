function psog_initialguess(beta0::T, dt0::T, Tfinal0::T, tol::T) where{T}
    # Derived parameters
    beta = beta0 / T(2)
    dt = dt0^2
    Tfinal = Tfinal0^2
    reps = tol
    delta = dt / Tfinal

    # Compute h with arbitrary precision
    h = T(1.5) * 2 * T(π) / (log(T(3)) + abs(beta) * log(1 / cos(T(1))) + log(1 / reps))

    # Compute integration limits
    tlower = (1 / abs(beta)) * log(reps * gamma(1 + beta))

    if beta >= 1
        tupper = log(1 / delta) + log(log(1 / reps)) + log(beta) + T(0.5)
    else
        tupper = log(1 / delta) + log(log(1 / reps))
    end

    # Integer bounds for summation
    M = floor(Int, tlower / h)
    N = ceil(Int, tupper / h)

    n1 = M:-1
    xs1 = -exp.(h .* T.(n1))
    ws1 = h / gamma(beta) .* exp.(beta * h .* T.(n1))
    @info "M, N: $M, $N" typeof(xs1) typeof(ws1)

    ws1new, xs1new = prony(xs1, ws1, tol)

    # matlab results
    # ws1new = [0.1297536498684488, 0.1138109239757789, 0.09983078384815716, 0.08781612517532403, 0.08102695587784227, 0.08709008482967978, 0.1004825231828261, 0.1121542677969847, 0.1200562218219118, 0.1240084595101501]
    # xs1new = [-0.7693582923088651, -0.5919121795202147, -0.4553916881934084, -0.35029070128222, -0.2680487783257383, -0.1973409274581753, -0.1304655035382046, -0.07114066567626039, -0.02674408299599483, -0.003035530072314068]

    @info length(xs1new)
    n2 = 0:N
    xs2 = -exp.(h .* T.(n2))
    ws2 = h / gamma(beta) .* exp.(beta * h .* T.(n2))

    xs = vcat(-real.(xs1new), -real.(xs2))  # vertical concatenation
    ws = vcat(real.(ws1new), real.(ws2))

    xs = xs ./ Tfinal
    ws = ws ./ Tfinal^beta

    nexp = length(ws)
    @info "Before model reduction, nexp = $(nexp)"

    return xs, ws
end

