@testset "DOFs.jl — Integrales de DOF tipo cara" begin
    mono = NedelecHex.Monomials.build_h2_monomials()

    # ── Consistencia con las filas 25–48 de G ───────────────────────────────
    _, G = NedelecHex.CoefficientGenerator.get_coefficients()
    face_block = zeros(24, 54)
    for j in 1:54, f in 1:6
        base = 24 + 4 * (f - 1)
        face_block[base-23:base-20, j] = NedelecHex.DOFs.evaluate_face_dof(mono[j], f)
    end
    @test isapprox(face_block, G[25:48, :], atol=1e-12)

    # ── Incompatibilidad direccional ────────────────────────────────────────
    # La cara 1 (z=0) sólo usa vectores de prueba Q_FACE en el plano xy (sin
    # componente ζ̂, ver DOFs.jl). Un campo puramente ζ̂ da los 4 DOFs nulos
    # sin importar su dependencia funcional.
    probe_zeta = p -> [0.0, 0.0, 1.0 + p[1] + p[2]^2 * p[3]]
    @test isapprox(NedelecHex.DOFs.evaluate_face_dof(probe_zeta, 1), zeros(4), atol=1e-12)

    # ── Coherencia de la suma de filas ──────────────────────────────────────
    # Q_FACE[1][1] + Q_FACE[1][2] = (1−x) + x = 1  ⇒  g₁(u)+g₂(u) = ∫_face u_x dA
    # Q_FACE[1][3] + Q_FACE[1][4] = (1−y) + y = 1  ⇒  g₃(u)+g₄(u) = ∫_face u_y dA
    probe = p -> [p[1]^2 + p[2], p[1] * p[2] + p[3], 0.0]
    dofs  = NedelecHex.DOFs.evaluate_face_dof(probe, 1)

    w2, p2 = NedelecHex.GaussQuadrature.gauss_quadrilateral(4)
    int_ux = sum(w2[k] * probe(NedelecHex.CoordinateTransforms.transform_2d_to_3d(p2[k,:], 1))[1]
                 for k in eachindex(w2))
    int_uy = sum(w2[k] * probe(NedelecHex.CoordinateTransforms.transform_2d_to_3d(p2[k,:], 1))[2]
                 for k in eachindex(w2))

    @test isapprox(dofs[1] + dofs[2], int_ux, atol=1e-12)
    @test isapprox(dofs[3] + dofs[4], int_uy, atol=1e-12)
end
