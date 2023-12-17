@testset "VP kernel accuracy, exp(-x^2)" begin
    f = x -> exp(-x^2)

    s, w = VP_cal(f, 4.0, 30, 100)
    
    @test max_error(f, s, w) < 1e-8
end

@testset "VP kernel accuracy, erfc(x)" begin
    f = x -> erfc(x)

    s, w = VP_cal(f, 4.0, 30, 100)

    @test max_error(f, s, w) < 1e-8
end