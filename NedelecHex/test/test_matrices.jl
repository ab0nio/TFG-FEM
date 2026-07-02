@testset "Assembly.jl — Matrices de rigidez y masa" begin
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    mesh    = NedelecHex.Mesh.box_mesh(nx=2, ny=2, nz=2)
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
    K, M    = NedelecHex.Assembly.assemble_global_matrices(
                  mesh, dof_map, W, mono, curl_mono)

    # ── Conteo analítico de DOFs para una malla estructurada nx=ny=nz=2 ──────
    # Nodos: (nx+1)(ny+1)(nz+1) = 3·3·3 = 27
    # Aristas de la rejilla (por eje): nx(ny+1)(nz+1) + (nx+1)ny(nz+1) + (nx+1)(ny+1)nz
    #                                = 2·3·3·3 = 54 aristas → 108 DOFs de arista
    # Caras de la rejilla (por eje):  (nx+1)ny·nz + nx(ny+1)nz + nx·ny(nz+1)
    #                                = 3·12       = 36 caras  → 144 DOFs de cara
    # Volumen: 8 elementos × 6 DOFs = 48 DOFs de volumen
    # Total: 108 + 144 + 48 = 300
    nx = ny = nz = 2
    n_nodes_expected = (nx+1) * (ny+1) * (nz+1)
    n_edges_expected = nx*(ny+1)*(nz+1) + (nx+1)*ny*(nz+1) + (nx+1)*(ny+1)*nz
    n_faces_expected = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1)
    n_elem_expected  = nx * ny * nz

    @test size(mesh.nodes, 1) == n_nodes_expected
    @test size(mesh.elements, 1) == n_elem_expected
    @test dof_map.n_dofs == 2*n_edges_expected + 4*n_faces_expected + 6*n_elem_expected
    @test NedelecHex.GlobalDOFs.n_edge_dofs(dof_map) == 2 * n_edges_expected
    @test NedelecHex.GlobalDOFs.n_face_dofs(dof_map) == 4 * n_faces_expected
    @test NedelecHex.GlobalDOFs.n_vol_dofs(dof_map)  == 6 * n_elem_expected

    dir  = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
    free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)

    # ── Densidad dispersa coherente con la conectividad esperada ─────────────
    # Cada elemento aporta como máximo 54×54 entradas COO; tras sumar
    # duplicados en los DOFs compartidos, nnz(K) no puede superar ese total.
    @test nnz(K) <= size(mesh.elements, 1) * 54 * 54
    @test nnz(K) == nnz(M)

    # ── Simetría ─────────────────────────────────────────────────────────────
    @test norm(K - K') / norm(K) < 1e-12
    @test norm(M - M') / norm(M) < 1e-12

    # ── Masa: definida positiva en DOFs libres ────────────────────────────────
    λM = eigvals(Symmetric(Matrix(M[free, free])))
    @test minimum(λM) > 0

    # ── Rigidez: semi-definida positiva en DOFs libres ────────────────────────
    # La energía de curl es no negativa: ∫|∇×u|² ≥ 0. El espacio nulo de K
    # corresponde a los modos irrotacionales (gradiente), que son espurios.
    λK = eigvals(Symmetric(Matrix(K[free, free])))
    @test minimum(λK) > -1e-8

    # ── Norma de la matriz elemental para el hexaedro unitario ───────────────
    # Valor de referencia calculado analíticamente con cuadratura de orden 4.
    nodes1 = [NedelecHex.ReferenceHex.NODES[k] for k in 1:8]
    Ke = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes1)
    Me = NedelecHex.Assembly.mass_matrix_element(W, mono, nodes1)
    @test isapprox(norm(Ke), 556.045466, rtol=1e-5)
    @test isapprox(norm(Me),  16.626423, rtol=1e-5)
end
