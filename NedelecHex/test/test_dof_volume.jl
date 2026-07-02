@testset "DOFs.jl — Integrales de DOF tipo volumen" begin
    mono = NedelecHex.Monomials.build_h2_monomials()

    # ── Consistencia con las filas 49–54 de G ───────────────────────────────
    _, G = NedelecHex.CoefficientGenerator.get_coefficients()
    vol_block = zeros(6, 54)
    for j in 1:54
        vol_block[:, j] = NedelecHex.DOFs.evaluate_volume_dof(mono[j])
    end
    @test isapprox(vol_block, G[49:54, :], atol=1e-12)

    # ── Valor analítico directo ──────────────────────────────────────────────
    # u = (1,0,0) constante sobre el hexaedro unidad:
    #   g_H,1(u) = ∫∫∫ (1−x) dV = 1/2 ,  g_H,2(u) = ∫∫∫ x dV = 1/2
    #   g_H,3..6(u) = 0  (u no tiene componente y ni z)
    const_field = p -> [1.0, 0.0, 0.0]
    dofs = NedelecHex.DOFs.evaluate_volume_dof(const_field)
    @test isapprox(dofs, [0.5, 0.5, 0.0, 0.0, 0.0, 0.0], atol=1e-12)

    # Mismo razonamiento para componentes η̂ y ζ̂
    dofs_y = NedelecHex.DOFs.evaluate_volume_dof(p -> [0.0, 1.0, 0.0])
    @test isapprox(dofs_y, [0.0, 0.0, 0.5, 0.5, 0.0, 0.0], atol=1e-12)

    dofs_z = NedelecHex.DOFs.evaluate_volume_dof(p -> [0.0, 0.0, 1.0])
    @test isapprox(dofs_z, [0.0, 0.0, 0.0, 0.0, 0.5, 0.5], atol=1e-12)
end
