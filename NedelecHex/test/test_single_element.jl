@testset "Assembly.jl + GlobalDOFs.jl — Ensamblaje de un solo elemento" begin
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    mesh    = NedelecHex.Mesh.box_mesh(nx=1, ny=1, nz=1)
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)

    @test size(mesh.elements, 1) == 1
    @test dof_map.n_dofs == 54

    K, M = NedelecHex.Assembly.assemble_global_matrices(
               mesh, dof_map, W, mono, curl_mono)

    @test size(K) == (54, 54)
    @test size(M) == (54, 54)

    # ── Simetría ─────────────────────────────────────────────────────────────
    @test norm(Matrix(K) - Matrix(K)') < 1e-12
    @test norm(Matrix(M) - Matrix(M)') < 1e-12

    # ── M definida positiva, K semidefinida positiva ─────────────────────────
    λM = eigvals(Symmetric(Matrix(M)))
    @test minimum(λM) > 0

    λK = eigvals(Symmetric(Matrix(K)))
    @test minimum(λK) > -1e-8

    # ── Dimensión del núcleo de K (modos gradiente espurios) ─────────────────
    # El núcleo de K en un único elemento corresponde a los campos
    # irrotacionales representables en H₂: los gradientes de la base nodal
    # triquadrática H¹ asociada (8 vértices + 12 puntos medios de arista +
    # 6 centros de cara + 1 centro de volumen = 27 nodos), menos la función
    # constante (gradiente nulo) → 27 − 1 = 26 modos.
    n_null = count(x -> abs(x) < 1e-6, λK)
    @test n_null == 26
end
