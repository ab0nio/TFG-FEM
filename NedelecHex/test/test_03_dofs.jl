@testset "DOFs.jl — Integrales de DOF tipo arista" begin
    # ── Integral analítica con campo no trivial ─────────────────────────────
    # Arista 1: de (0,0,0) a (1,0,0), tangente ξ̂, parámetro t = ξ, longitud 1.
    # Sobre la arista η = ζ = 0, así que u·τ̂ = u_ξ = ξ² = t².
    #   g₁(u) = ∫₀¹ t²(1−t) dt = 1/3 − 1/4 = 1/12
    #   g₂(u) = ∫₀¹ t²·t    dt = 1/4
    probe_quad = p -> [p[1]^2 + p[2], p[3], p[1]]
    dofs_quad = NedelecHex.DOFs.evaluate_edge_dof(probe_quad, 1)
    explica("Arista 1 con u_ξ = t² sobre la arista: a mano ∫t²(1−t)dt = 1/12 y ∫t³dt = 1/4; el código da $(dofs_quad). Al ser valores distintos entre sí, un intercambio de q₁=1−t con q₂=t no pasaría desapercibido.")
    @test isapprox(dofs_quad, [1/12, 1/4], atol=1e-12)

    # ── Ortogonalidad direccional ────────────────────────────────────────────
    # Las aristas 1 y 3 (bordes inferior/superior en x) son tangentes a ξ̂.
    # Un campo con componente únicamente en η̂ o ζ̂ es siempre perpendicular a
    # esa tangente, así que su DOF de arista debe ser exactamente cero.
    explica("Ortogonalidad direccional: las aristas 1 y 3 son tangentes a ξ̂; un campo puramente η̂ o ζ̂ tiene componente tangencial nula sobre ellas ⇒ DOF = [0,0] exacto.")
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
    explica("Anulación geométrica: la arista 9 vive en ξ=η=0; el campo u_ζ = ξ se anula sobre toda ella ⇒ DOF nulo aunque el campo apunte en la dirección tangente ζ̂.")
    probe_vanish = p -> [0.0, 0.0, p[1]]
    @test isapprox(NedelecHex.DOFs.evaluate_edge_dof(probe_vanish, 9), [0.0, 0.0], atol=1e-12)

    # ── Valor analítico directo ──────────────────────────────────────────────
    # u = (1,0,0) constante sobre la arista 1 (ξ̂, longitud 1):
    #   g₁(u) = ∫₀¹ 1·(1−t) dt = 1/2 ,  g₂(u) = ∫₀¹ 1·t dt = 1/2
    explica("Campo constante u = (1,0,0) sobre la arista 1: ∫(1−t)dt = ∫t dt = 1/2 ⇒ DOFs = [0.5, 0.5].")
    const_field = p -> [1.0, 0.0, 0.0]
    @test isapprox(NedelecHex.DOFs.evaluate_edge_dof(const_field, 1), [0.5, 0.5], atol=1e-12)
end
@testset "DOFs.jl — Integrales de DOF tipo cara" begin
    # ── Integral analítica con campo no trivial ─────────────────────────────
    # Cara 1 (z=0): g_qi(u) = ∫₀¹∫₀¹ u·q dxdy con q ∈ {(1−x)x̂, x·x̂, (1−y)ŷ, y·ŷ}.
    # Con u = (xy, x + y², ·) — la componente ζ̂ no interviene:
    #   g₁ = ∫∫ xy(1−x)      = (1/2−1/3)·(1/2)      = 1/12
    #   g₂ = ∫∫ x²y          = (1/3)·(1/2)          = 1/6
    #   g₃ = ∫∫ (x+y²)(1−y)  = (1/2)(1/2) + (1/3−1/4) = 1/3
    #   g₄ = ∫∫ (x+y²)y      = (1/2)(1/2) + 1/4       = 1/2
    probe_face = p -> [p[1] * p[2], p[1] + p[2]^2, 7.0 * p[3] + 1.0]
    dofs_face = NedelecHex.DOFs.evaluate_face_dof(probe_face, 1)
    explica("Cara 1 (z=0) con u = (xy, x+y², ·): integrales a mano [1/12, 1/6, 1/3, 1/2]; el código da $(dofs_face). El campo lleva componente ζ̂ a propósito para verificar que la cara la ignora.")
    @test isapprox(dofs_face, [1/12, 1/6, 1/3, 1/2], atol=1e-12)

    # ── Incompatibilidad direccional ────────────────────────────────────────
    # La cara 1 (z=0) sólo usa vectores de prueba Q_FACE en el plano xy (sin
    # componente ζ̂, ver DOFs.jl). Un campo puramente ζ̂ da los 4 DOFs nulos
    # sin importar su dependencia funcional.
    explica("Incompatibilidad direccional: los vectores de prueba de la cara z=0 viven en el plano xy ⇒ un campo puramente ζ̂ da los 4 DOFs nulos sea cual sea su dependencia funcional.")
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

    explica("Coherencia de sumas: (1−x)+x = 1 ⇒ g₁+g₂ debe igualar ∫u_x dA = $(int_ux) (suma de DOFs: $(dofs[1]+dofs[2])) y g₃+g₄ debe igualar ∫u_y dA = $(int_uy) (suma: $(dofs[3]+dofs[4])).")
    @test isapprox(dofs[1] + dofs[2], int_ux, atol=1e-12)
    @test isapprox(dofs[3] + dofs[4], int_uy, atol=1e-12)
end
@testset "DOFs.jl — Integrales de DOF tipo volumen" begin
    # ── Integral analítica con campo no trivial ─────────────────────────────
    # g_qi(u) = ∫∫∫ u·q dV sobre [0,1]³, q ∈ {(1−x)x̂, x·x̂, (1−y)ŷ, y·ŷ, (1−z)ẑ, z·ẑ}.
    # Con u = (xy, z, x²), integrales separables (∫x=1/2, ∫x²=1/3, ∫x(1−x)=1/6):
    #   g₁ = ∫∫∫ xy(1−x) = (1/6)(1/2) = 1/12    g₂ = ∫∫∫ x²y    = (1/3)(1/2) = 1/6
    #   g₃ = ∫∫∫ z(1−y)  = (1/2)(1/2) = 1/4     g₄ = ∫∫∫ zy     = (1/2)(1/2) = 1/4
    #   g₅ = ∫∫∫ x²(1−z) = (1/3)(1/2) = 1/6     g₆ = ∫∫∫ x²z    = (1/3)(1/2) = 1/6
    probe_vol = p -> [p[1] * p[2], p[3], p[1]^2]
    dofs_vol = NedelecHex.DOFs.evaluate_volume_dof(probe_vol)
    explica("Volumen con u = (xy, z, x²): integrales separables a mano [1/12, 1/6, 1/4, 1/4, 1/6, 1/6]; el código da $(dofs_vol).")
    @test isapprox(dofs_vol, [1/12, 1/6, 1/4, 1/4, 1/6, 1/6], atol=1e-12)

    # ── Valor analítico directo ──────────────────────────────────────────────
    # u = (1,0,0) constante sobre el hexaedro unidad:
    #   g_H,1(u) = ∫∫∫ (1−x) dV = 1/2 ,  g_H,2(u) = ∫∫∫ x dV = 1/2
    #   g_H,3..6(u) = 0  (u no tiene componente y ni z)
    explica("Campos constantes unitarios: cada pareja Lagrange debe dar [0.5, 0.5] en su eje y 0 en los demás — verifica el orden (x,x,y,y,z,z) de Q_VOL.")
    const_field = p -> [1.0, 0.0, 0.0]
    dofs = NedelecHex.DOFs.evaluate_volume_dof(const_field)
    @test isapprox(dofs, [0.5, 0.5, 0.0, 0.0, 0.0, 0.0], atol=1e-12)

    # Mismo razonamiento para componentes η̂ y ζ̂
    dofs_y = NedelecHex.DOFs.evaluate_volume_dof(p -> [0.0, 1.0, 0.0])
    @test isapprox(dofs_y, [0.0, 0.0, 0.5, 0.5, 0.0, 0.0], atol=1e-12)

    dofs_z = NedelecHex.DOFs.evaluate_volume_dof(p -> [0.0, 0.0, 1.0])
    @test isapprox(dofs_z, [0.0, 0.0, 0.0, 0.0, 0.5, 0.5], atol=1e-12)
end
