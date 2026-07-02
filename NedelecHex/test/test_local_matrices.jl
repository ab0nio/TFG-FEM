@testset "Assembly.jl — Matrices locales Kᵉ, Mᵉ en el elemento de referencia" begin
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()
    nodes     = NedelecHex.ReferenceHex.NODES

    Ke = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes)
    Me = NedelecHex.Assembly.mass_matrix_element(W, mono, nodes)

    @test size(Ke) == (54, 54)
    @test size(Me) == (54, 54)

    # ── Simetría ─────────────────────────────────────────────────────────────
    @test norm(Ke - Ke') < 1e-12
    @test norm(Me - Me') < 1e-12

    # ── Mᵉ definida positiva ─────────────────────────────────────────────────
    λM = eigvals(Symmetric(Me))
    @test minimum(λM) > 0

    # ── Kᵉ semidefinida positiva ──────────────────────────────────────────────
    λK = eigvals(Symmetric(Ke))
    @test minimum(λK) > -1e-8

    # ── Número de condicionamiento (Tabla III del paper: 19 y 35) ────────────
    # El paper reporta el condicionamiento de las matrices *precondicionadas*
    # con un escalado diagonal de Jacobi (Mᵖ = D⁻¹ᐟ²MD⁻¹ᐟ², análogo para Kᵖ).
    # Para Kᵖ el precondicionador sólo tiene sentido sobre el complemento del
    # espacio nulo (autovalores espurios ≈0 se excluyen).
    Dm = Diagonal(1.0 ./ sqrt.(diag(Me)))
    Mp = Dm * Me * Dm
    @test isapprox(cond(Mp), 19.0, atol=0.5)

    Dk       = Diagonal(1.0 ./ sqrt.(diag(Ke) .+ (diag(Ke) .== 0)))
    Kp       = Dk * Ke * Dk
    λKp_nz   = filter(x -> x > 1e-6, eigvals(Symmetric(Kp)))
    @test isapprox(maximum(λKp_nz) / minimum(λKp_nz), 35.0, atol=1.0)
end
