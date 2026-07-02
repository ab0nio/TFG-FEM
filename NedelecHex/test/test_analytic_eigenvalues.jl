@testset "HalfFilledCavity.jl — Eigenvalores analíticos de la cavidad semillena" begin
    kc = NedelecHex.HalfFilledCavity.analytic_eigenvalues()

    # Parámetros del benchmark (paper Amor-Martin 2023):
    #   cavidad  a=1cm, b=1mm, c=1cm  |  losa dieléctrica h=5mm, εᵣ=2
    # Los 6 modos TE₁₀ₙ / TE₂₀ₙ más bajos deben estar ordenados y coincidir
    # con los valores tabulados en el paper con tolerancia 0.5 rad/m.
    @test length(kc) == 6
    @test all(kc .> 0)
    @test issorted(kc)

    expected = [353.78, 544.50, 599.80, 750.31, 763.32, 809.58]
    labels   = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        @test isapprox(kc[i], expected[i], atol=0.5)
    end

    for (label, k) in zip(labels, kc)
        @printf("%-6s = %10.4f rad/m\n", label, k)
    end

    # ── Convergencia de la bisección ─────────────────────────────────────────
    # Cada kc encontrada debe ser una raíz genuina de la ecuación trascendente:
    # el residuo debe ser prácticamente nulo y debe existir un cambio de signo
    # real a ambos lados (no un artefacto cerca de una asíntota de cot()).
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    modes = [(1,0), (2,0), (1,0), (3,0), (2,0), (1,0)]
    for (i, k) in enumerate(kc)
        m, n = modes[i]
        residual = NedelecHex.HalfFilledCavity._resonance_eq(k, m, n, a, b, c, h, epr)
        @test abs(residual) < 1e-6

        r_lo = NedelecHex.HalfFilledCavity._resonance_eq(k - 1e-6, m, n, a, b, c, h, epr)
        r_hi = NedelecHex.HalfFilledCavity._resonance_eq(k + 1e-6, m, n, a, b, c, h, epr)
        @test sign(r_lo) != sign(r_hi)
    end
end
