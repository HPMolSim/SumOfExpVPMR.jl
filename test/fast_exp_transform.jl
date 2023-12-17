@testset "test fast exponentials transform 1d" begin
    f = x -> exp(-x^2)
    soepara = VPMR_cal(f, 6.0, 60, 200, 16)[1]
    for T in [Float64, ComplexF64]
        q_1 = rand(T, 100)
        q_2 = rand(T, 100)
        x = rand(100)
        α = rand()

        sum_exact = zero(ComplexF64)
        for i in 1:length(x)
            for j in 1:length(x)
                sum_exact += q_1[i] * q_2[j] * f(abs(x[i] - x[j])/ α)
            end
        end

        sum_soe = FET1d(q_1, q_2, x, soepara, α = α)

        @test isapprox((sum_exact), (sum_soe), rtol = 1e-8)
    end
end