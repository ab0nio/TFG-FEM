@testset "DOFs.jl: Integrales de DOF tipo arista" begin
    # Arista 1: de (0,0,0) a (1,0,0), tangente ξ̂, parámetro t = ξ.
    # Sobre la arista η=ζ=0, así que u·τ̂ = u_ξ = ξ² = t².
    # g₁(u) = ∫₀¹ t²(1−t) dt = 1/12,   g₂(u) = ∫₀¹ t³ dt = 1/4
    probe = p -> [p[1]^2 + p[2], p[3], p[1]]
    dofs = NedelecHex.DOFs.evaluate_edge_dof(probe, 1)
    explica("Arista 1 con u_ξ = t² sobre la arista: a mano ∫t²(1−t)dt = 1/12 y ∫t³dt = 1/4; el código da $(dofs).")
    @test isapprox(dofs, [1/12, 1/4], atol=1e-12)

    # Las aristas 1 y 3 son tangentes a ξ̂: un campo sin componente ξ̂ da DOF
    # nulo por ortogonalidad. La arista 9 (ξ=η=0) anula un campo con
    # componente ζ̂ que se anula geométricamente sobre ella (u_ζ = ξ).
    explica("COmprueba Dos formas de obtener un DOF nulo: por ortogonalidad ,en las aristas 1 y 3 ,
    por ejemplo (dirección ξ̂) un campo sin componente ξ̂ se proyecta a cero
    ó por anulación geométrica, por ejemplo ,la arista 9 tiene ξ=η=0, y ahí el campo u_ζ=ξ,
    que no es nulo en el cubo en general, vale cero justo sobre esa arista.
    Un test comprueba que integra bien y el otro que estamos integrando en la arista correcta")
    probe_perp = p -> [0.0, p[1] + p[2]^2, p[3]^2]
    for e in (1, 3)
        @test NedelecHex.DOFs.evaluate_edge_dof(probe_perp, e) == [0.0, 0.0]
    end
    probe_vanish = p -> [0.0, 0.0, p[1]]
    @test isapprox(NedelecHex.DOFs.evaluate_edge_dof(probe_vanish, 9), [0.0, 0.0], atol=1e-12)
end

@testset "DOFs.jl: Integrales de DOF tipo cara" begin
    # Cara 1 (z=0): g_qi(u) = ∫₀¹∫₀¹ u·q dxdy con q ∈ {(1−x)x̂, x·x̂, (1−y)ŷ, y·ŷ}.
    # u = (xy, x+y², 7z+1): la componente ζ̂ se incluye a propósito para
    # comprobar que la cara la ignora.
    # g₁=1/12, g₂=1/6, g₃=1/3, g₄=1/2
    probe = p -> [p[1] * p[2], p[1] + p[2]^2, 7.0 * p[3] + 1.0]
    dofs = NedelecHex.DOFs.evaluate_face_dof(probe, 1)
    explica("Comprueba que en Cara 1 (z=0) con u=(xy, x+y², 7z+1) debe salir:[1/12,1/6,1/3,1/2] y el código da: $(dofs).")
    @test isapprox(dofs, [1/12, 1/6, 1/3, 1/2], atol=1e-12)

    # Campo puramente ζ̂ ⇒ los 4 DOFs de la cara z=0 (vectores de prueba en el
    # plano xy) son nulos sea cual sea su dependencia funcional.
    probe_zeta = p -> [0.0, 0.0, 1.0 + p[1] + p[2]^2 * p[3]]
    @test isapprox(NedelecHex.DOFs.evaluate_face_dof(probe_zeta, 1), zeros(4), atol=1e-12)

    # Q_FACE[1][1]+Q_FACE[1][2] = (1−x)+x = 1 y Q_FACE[1][3]+Q_FACE[1][4] = (1−y)+y = 1,
    # así que g₁+g₂ = ∫u_x dA y g₃+g₄ = ∫u_y dA. Con u = (x²+y, xy, 0):
    # ∫∫(x²+y) = 1/3+1/2 = 5/6,   ∫∫xy = 1/4
    probe2 = p -> [p[1]^2 + p[2], p[1] * p[2], 0.0]
    dofs2 = NedelecHex.DOFs.evaluate_face_dof(probe2, 1)
    explica("La suma del DOF 1 con DOF2 debe ser igual a la integral de la componente en X(lo mismo con Y)
    ya que :Q_FACE[1][1]+Q_FACE[1][2] = (1−x)+x = 1 y Q_FACE[1][3]+Q_FACE[1][4] = (1−y)+y = 1, entonces
    g₁+g₂ debe ser ∫u_x dA=5/6 y g₃+g₄=∫u_y dA=1/4; código da $(dofs2[1]+dofs2[2]) y $(dofs2[3]+dofs2[4]).")
    @test isapprox(dofs2[1] + dofs2[2], 5/6, atol=1e-12)
    @test isapprox(dofs2[3] + dofs2[4], 1/4, atol=1e-12)
end

@testset "DOFs.jl: Integrales de DOF tipo volumen" begin
    # g_qi(u) = ∫∫∫ u·q dV sobre [0,1]³, q ∈ {(1−x)x̂, x·x̂, (1−y)ŷ, y·ŷ, (1−z)ẑ, z·ẑ}.
    # u = (xy, z, x²): g₁=1/12, g₂=1/6, g₃=1/4, g₄=1/4, g₅=1/6, g₆=1/6
    probe = p -> [p[1] * p[2], p[3], p[1]^2]
    dofs = NedelecHex.DOFs.evaluate_volume_dof(probe)
    explica("Comprueba que en Volumen con u=(xy,z,x²): debe dar [1/12,1/6,1/4,1/4,1/6,1/6] y en el código da $(dofs).")
    @test isapprox(dofs, [1/12, 1/6, 1/4, 1/4, 1/6, 1/6], atol=1e-12)

    # Campo constante unitario en cada eje: debe dar [0.5,0.5] en su pareja
    # Lagrange y 0 en las demás: verifica el orden (x,x,y,y,z,z) de Q_VOL.
    explica("Comprueba que cada uno da [0.5,0.5] en su pareja y 0 en las demás.")
    for (axis, expected) in ((1, [0.5, 0.5, 0.0, 0.0, 0.0, 0.0]),
        (2, [0.0, 0.0, 0.5, 0.5, 0.0, 0.0]),
        (3, [0.0, 0.0, 0.0, 0.0, 0.5, 0.5]))
        field = p -> [axis == 1 ? 1.0 : 0.0, axis == 2 ? 1.0 : 0.0, axis == 3 ? 1.0 : 0.0]
        @test isapprox(NedelecHex.DOFs.evaluate_volume_dof(field), expected, atol=1e-12)
    end
end
