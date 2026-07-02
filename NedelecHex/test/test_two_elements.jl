@testset "GlobalDOFs.jl + Assembly.jl — Ensamblaje de dos elementos adyacentes" begin
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    mesh    = NedelecHex.Mesh.box_mesh(nx=2, ny=1, nz=1)
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)

    # ── Conteo de DOFs ────────────────────────────────────────────────────────
    # Dos elementos de 54 DOFs cada uno comparten la cara x=0.5:
    #   4 aristas × 2 DOFs = 8 DOFs de arista compartidos
    #   1 cara    × 4 DOFs = 4 DOFs de cara   compartidos
    #   Total compartidos  = 12  →  96 DOFs globales (54+54−12)
    @test dof_map.n_dofs == 96

    conn1  = dof_map.connectivity[1, :]
    conn2  = dof_map.connectivity[2, :]
    shared = sort!(collect(intersect(Set(conn1), Set(conn2))))
    @test length(shared) == 12

    shared_edge = intersect(Set(conn1[1:24]),  Set(conn2[1:24]))
    shared_face = intersect(Set(conn1[25:48]), Set(conn2[25:48]))
    shared_vol  = intersect(Set(conn1[49:54]), Set(conn2[49:54]))
    @test length(shared_edge) == 8    # 4 aristas × 2 DOFs
    @test length(shared_face) == 4    # 1 cara × 4 DOFs
    @test isempty(shared_vol)

    # ── Acumulación en DOFs compartidos ──────────────────────────────────────
    # El ensamblaje COO suma automáticamente las contribuciones de ambos
    # elementos cuando tienen el mismo índice global. Para cada DOF compartido g:
    #   K_global[g, g]  =  K1[loc1, loc1]  +  K2[loc2, loc2]
    nidx1  = mesh.elements[1, :]
    nidx2  = mesh.elements[2, :]
    nodes1 = [mesh.nodes[nidx1[k], :] for k in 1:8]
    nodes2 = [mesh.nodes[nidx2[k], :] for k in 1:8]

    K1 = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes1)
    K2 = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes2)
    K, _ = NedelecHex.Assembly.assemble_global_matrices(
               mesh, dof_map, W, mono, curl_mono)

    for g in shared
        loc1 = findfirst(==(g), conn1)
        loc2 = findfirst(==(g), conn2)
        k1v  = isnothing(loc1) ? 0.0 : K1[loc1, loc1]
        k2v  = isnothing(loc2) ? 0.0 : K2[loc2, loc2]
        @test isapprox(K[g, g], k1v + k2v, atol=1e-10)
    end

    # Lo mismo pero para entradas fuera de la diagonal entre dos DOFs
    # compartidos cualesquiera: K_global[g,h] = K1[loc1g,loc1h] + K2[loc2g,loc2h]
    # (0 si el DOF local no existe en ese elemento).
    for g in shared, h in shared
        g == h && continue
        loc1g, loc1h = findfirst(==(g), conn1), findfirst(==(h), conn1)
        loc2g, loc2h = findfirst(==(g), conn2), findfirst(==(h), conn2)
        k1v = (isnothing(loc1g) || isnothing(loc1h)) ? 0.0 : K1[loc1g, loc1h]
        k2v = (isnothing(loc2g) || isnothing(loc2h)) ? 0.0 : K2[loc2g, loc2h]
        @test isapprox(K[g, h], k1v + k2v, atol=1e-10)
    end

    # ── Orientación consistente de las 4 aristas de la cara compartida ───────
    # build_dof_map asigna, para cada arista física, un DOF "nodo mínimo" y un
    # DOF "nodo máximo" (canónicos, según minmax(nodo_global)), permutando los
    # dos DOFs locales cuando la arista local corre en sentido contrario. Aquí
    # se recupera de forma independiente, para cada arista compartida, qué DOF
    # global quedó ligado al nodo mínimo y al máximo según CADA elemento por
    # separado, y se verifica que ambos elementos coinciden — es decir, que la
    # orientación relativa (ξ̂ o η̂ o ζ̂, según corresponda) es consistente a
    # través de la cara compartida.
    EDGES = NedelecHex.ReferenceHex.EDGES
    n_shared_edges_checked = 0
    for (le1, (ln1, ln2)) in enumerate(EDGES)
        gn1, gn2 = nidx1[ln1], nidx1[ln2]
        for (le2, (lm1, lm2)) in enumerate(EDGES)
            gm1, gm2 = nidx2[lm1], nidx2[lm2]
            Set([gn1, gn2]) != Set([gm1, gm2]) && continue

            base1 = 2 * (le1 - 1)
            base2 = 2 * (le2 - 1)
            dof_min1, dof_max1 = gn1 < gn2 ? (conn1[base1+1], conn1[base1+2]) :
                                              (conn1[base1+2], conn1[base1+1])
            dof_min2, dof_max2 = gm1 < gm2 ? (conn2[base2+1], conn2[base2+2]) :
                                              (conn2[base2+2], conn2[base2+1])

            @test dof_min1 == dof_min2
            @test dof_max1 == dof_max2
            n_shared_edges_checked += 1
        end
    end
    @test n_shared_edges_checked == 4

    # ── Matrices globales K y M: simetría y M definida positiva ─────────────
    _, M = NedelecHex.Assembly.assemble_global_matrices(
               mesh, dof_map, W, mono, curl_mono)
    @test norm(Matrix(K) - Matrix(K)') < 1e-12
    @test norm(Matrix(M) - Matrix(M)') < 1e-12
    @test minimum(eigvals(Symmetric(Matrix(M)))) > 0

    # ── Formato COO: sin duplicados, índices válidos ─────────────────────────
    # sparse() de SparseArrays suma automáticamente las entradas con (i,j)
    # repetido, así que la propia existencia de K,M como SparseMatrixCSC bien
    # formadas ya garantiza ausencia de duplicados. Se comprueba además que
    # todos los índices de fila/columna están en rango.
    rows, cols, _ = findnz(K)
    @test all(1 .<= rows .<= dof_map.n_dofs)
    @test all(1 .<= cols .<= dof_map.n_dofs)
    @test length(rows) == length(unique(zip(rows, cols)))

    @info "Tabla de conectividad local→global (dos elementos adyacentes)" conn1 conn2
end
