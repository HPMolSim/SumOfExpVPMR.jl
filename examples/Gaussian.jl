# Approximating Gaussian

using Plots, SumOfExpVPMR

begin
    f = x -> exp(-x^2)
    
    # 4 terms approximation
    sw4, σ4 = VPMR_cal(f, 4.0, 40, 100, 4, print_info = true)
    error4 = soe_error(f, sw4)

    # 8 terms approximation
    sw8, σ8 = VPMR_cal(f, 4.0, 40, 100, 8, print_info = true)
    error8 = soe_error(f, sw8)

    # 12 terms approximation
    sw12, σ12 = VPMR_cal(f, 4.0, 40, 100, 12, print_info = true)
    error12 = soe_error(f, sw12)

    # 16 terms approximation
    sw16, σ16 = VPMR_cal(f, 6.0, 60, 200, 16, print_info = true)
    error16 = soe_error(f, sw16)
end

begin
    x = [0.0:0.01:10.0...]
    plot(xlabel = "x", ylabel = "log(error)", dpi = 500, title = "Approximating Gaussian via soe")
    plot!(x, log10.(error4), label = "M = 4")
    plot!(x, log10.(error8), label = "M = 8")
    plot!(x, log10.(error12), label = "M = 12")
    plot!(x, log10.(error16), label = "M = 16")
    savefig("Gaussian.png")
end