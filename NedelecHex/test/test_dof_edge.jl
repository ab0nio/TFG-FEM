@testset "DOFs.jl — Integrales de DOF tipo arista" begin
    mono = NedelecHex.Monomials.build_h2_monomials()

    # ── Consistencia con las filas 1–24 de G ────────────────────────────────
    # get_coefficients() construye G llamando exactamente a evaluate_edge_dof
    # por cada arista y monomio; reconstruirlo aquí de forma independiente y
    # compararlo confirma que DOFs.jl y CoefficientGenerator.jl están de
    # acuerdo.
    _, G = NedelecHex.CoefficientGenerator.get_coefficients()
    edge_block = zeros(24, 54)
    for j in 1:54, e in 1:12
        edge_block[2*e-1:2*e, j] = NedelecHex.DOFs.evaluate_edge_dof(mono[j], e)
    end
    @test isapprox(edge_block, G[1:24, :], atol=1e-12)

    # ── Ortogonalidad direccional ────────────────────────────────────────────
    # Las aristas 1 y 3 (bordes inferior/superior en x) son tangentes a ξ̂.
    # Un campo con componente únicamente en η̂ o ζ̂ es siempre perpendicular a
    # esa tangente, así que su DOF de arista debe ser exactamente cero.
    probe_eta  = p -> [0.0, p[1] + p[2]^2 + p[3], 0.0]
    probe_zeta = p -> [0.0, 0.0, p[1]^2 + p[2] + p[3]^2]
    for e in (1, 3)
        @test NedelecHex.DOFs.evaluate_edge_dof(probe_eta,  e) == [0.0, 0.0]
        @test NedelecHex.DOFs.evaluate_edge_dof(probe_zeta, e) == [0.0, 0.0]
    end

    # ── Anulación geométrica ─────────────────────────────────────────────────
    # La arista 9 va de (0,0,0) a (0,0,1): ξ=η=0 en todo punto de la arista.
    # Un campo ζ̂ que se anula cuando ξ=0 (p.ej. u_ζ = ξ) da un DOF nulo aunque
    # su componente esté alineada con la tangente de la arista (ζ̂).
    probe_vanish = p -> [0.0, 0.0, p[1]]
    @test isapprox(NedelecHex.DOFs.evaluate_edge_dof(probe_vanish, 9), [0.0, 0.0], atol=1e-12)

    # ── Valor analítico directo ──────────────────────────────────────────────
    # u = (1,0,0) constante sobre la arista 1 (ξ̂, longitud 1):
    #   g₁(u) = ∫₀¹ 1·(1−t) dt = 1/2 ,  g₂(u) = ∫₀¹ 1·t dt = 1/2
    const_field = p -> [1.0, 0.0, 0.0]
    @test isapprox(NedelecHex.DOFs.evaluate_edge_dof(const_field, 1), [0.5, 0.5], atol=1e-12)
end
