# Ensamblaje global, progresión 1 → 2 → 8 elementos. Cada nivel añade una
# dificultad sobre el anterior: un elemento (sin vecinos), dos elementos (primer
# reparto de DOFs compartidos, verificado al detalle) y una malla 2×2×2 (conteos
# combinatorios y propiedades globales). Prueba Assembly.jl y GlobalDOFs.jl.

# Malla estructurada + mapa de DOFs + matrices globales en un paso
# (los fixtures W_REF, MONO, CURL_MONO son globales de runtests.jl).
function assemble_box(nx, ny, nz)
    mesh = NedelecHex.Mesh.box_mesh(nx=nx, ny=ny, nz=nz)
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
    K, M = NedelecHex.Assembly.assemble_global_matrices(
        mesh, dof_map, W_REF, MONO, CURL_MONO)
    return mesh, dof_map, K, M
end

@testset "Assembly: 1 elemento" begin
    mesh, dof_map, K, M = assemble_box(1, 1, 1)

    @test size(mesh.elements, 1) == 1
    @test dof_map.n_dofs == 54
    explica("Malla de 1 elemento: los 54 DOFs locales son los 54 globales.")

    @test size(K) == (54, 54)
    @test size(M) == (54, 54)

    # Dimensión del núcleo de K (modos gradiente espurios)
    # El núcleo de K en un único elemento corresponde a los campos
    # irrotacionales representables en H₂: los gradientes de la base nodal
    # triquadrática H¹ asociada (8 vértices + 12 puntos medios de arista +
    # 6 centros de cara + 1 centro de volumen = 27 nodos), menos la función
    # constante (gradiente nulo) → 27 − 1 = 26 modos.
    λK = eigvals(Symmetric(Matrix(K)))
    n_null = count(x -> abs(x) < 1e-6, λK)
    explica("Dimensión del núcleo de K: $(n_null) autovalores ≈ 0
    Deben salir 26: (8 vértices + 12 puntos medios de arista +
    6 centros de cara + 1 centro de volumen = 27 nodos), menos la función
    constante da  27 − 1.")
    @test n_null == 26
end

@testset "Assembly: 2 elementos adyacentes (DOFs compartidos)" begin
    mesh, dof_map, K, _ = assemble_box(2, 1, 1)

    # Conteo de DOFs
    # Dos elementos de 54 DOFs comparten la cara x=0.5:
    # 4 aristas × 2 DOFs + 1 cara × 4 DOFs = 12 compartidos (los de volumen
    # nunca se comparten) → 96 DOFs globales (54 + 54 − 12).
    explica("Dos elementos que comparten la cara x=0.5: Deben salir :54+54−12 = 96 DOFs globales y hay: $(dof_map.n_dofs) ; los 12 compartidos son 4 aristas × 2 + 1 cara × 4 (los de volumen nunca se comparten).")
    @test dof_map.n_dofs == 96

    conn1 = dof_map.connectivity[1, :]
    conn2 = dof_map.connectivity[2, :]
    shared = sort!(collect(intersect(Set(conn1), Set(conn2))))
    @test length(shared) == 12
    @test length(intersect(Set(conn1[1:24]), Set(conn2[1:24]))) == 8   # 4 aristas × 2
    @test length(intersect(Set(conn1[25:48]), Set(conn2[25:48]))) == 4   # 1 cara × 4
    @test isempty(intersect(Set(conn1[49:54]), Set(conn2[49:54])))        # volumen: nunca

    # Acumulación COO en DOFs compartidos
    # El ensamblaje suma las contribuciones de ambos elementos cuando comparten
    # índice global: K_global[g,h] = K1[loc1] + K2[loc2] (0 si el DOF no está en
    # ese elemento). Se comprueba entrada a entrada; el caso g==h cubre la
    # diagonal.
    nidx1, nidx2 = mesh.elements[1, :], mesh.elements[2, :]
    nodes1, nodes2 = NedelecHex.Mesh.element_nodes(mesh, 1), NedelecHex.Mesh.element_nodes(mesh, 2)
    K1 = NedelecHex.Assembly.stiff_matrix_element(W_REF, CURL_MONO, nodes1)
    K2 = NedelecHex.Assembly.stiff_matrix_element(W_REF, CURL_MONO, nodes2)

    explica("para cada pareja de DOFs compartidos (g,h), K_global debe ser (k1 + k2).")
    for g in shared, h in shared
        loc1g, loc1h = findfirst(==(g), conn1), findfirst(==(h), conn1)
        loc2g, loc2h = findfirst(==(g), conn2), findfirst(==(h), conn2)
        k1 = (isnothing(loc1g) || isnothing(loc1h)) ? 0.0 : K1[loc1g, loc1h]
        k2 = (isnothing(loc2g) || isnothing(loc2h)) ? 0.0 : K2[loc2g, loc2h]
        @test isapprox(K[g, h], k1 + k2, atol=1e-10)
    end

    # Orientación consistente de las 4 aristas de la cara compartida
    # Para cada arista física compartida, build_dof_map liga un DOF al nodo
    # global mínimo y otro al máximo (permutando los locales si la arista corre
    # al revés). Ambos elementos deben coincidir en esa asignación: de ello
    # depende la continuidad tangencial a través de la cara.
    explica("Comprueba que en las 4 aristas compartidas la orientación es consistente(en las 4 aristas hay el mismo nodo mínimo y máximo).")
    EDGES = NedelecHex.ReferenceHex.EDGES
    n_checked = 0
    for (le1, (ln1, ln2)) in enumerate(EDGES), (le2, (lm1, lm2)) in enumerate(EDGES)
        gn1, gn2 = nidx1[ln1], nidx1[ln2]
        gm1, gm2 = nidx2[lm1], nidx2[lm2]
        Set([gn1, gn2]) != Set([gm1, gm2]) && continue

        b1, b2 = 2*(le1-1), 2*(le2-1)
        minmax1 = gn1 < gn2 ? (conn1[b1+1], conn1[b1+2]) : (conn1[b1+2], conn1[b1+1])
        minmax2 = gm1 < gm2 ? (conn2[b2+1], conn2[b2+2]) : (conn2[b2+2], conn2[b2+1])
        @test minmax1 == minmax2
        n_checked += 1
    end
    @test n_checked == 4

    TEST_VERBOSE && @info "Conectividad local→global (2 elementos adyacentes)" conn1 conn2
end

@testset "Assembly: malla 2×2×2" begin
    mesh, dof_map, K, M = assemble_box(2, 2, 2)

    # Conteo analítico de DOFs para una malla estructurada nx=ny=nz=2
    # Aristas: nx(ny+1)(nz+1)+(nx+1)ny(nz+1)+(nx+1)(ny+1)nz = 54 → 108 DOFs
    # Caras:   (nx+1)ny·nz + nx(ny+1)nz + nx·ny(nz+1)        = 36 → 144 DOFs
    # Volumen: 8 elementos × 6                               = 48  DOFs
    # Total: 108 + 144 + 48 = 300.
    nx = ny = nz = 2
    n_nodes = (nx+1)*(ny+1)*(nz+1)
    n_edges = nx*(ny+1)*(nz+1) + (nx+1)*ny*(nz+1) + (nx+1)*(ny+1)*nz
    n_faces = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1)
    n_elem = nx*ny*nz

    explica("Número de DOFs (malla 2×2×2): $(2*n_edges) DOFs de arista + $(4*n_faces) de cara + $(6*n_elem) de volumen = $(2*n_edges + 4*n_faces + 6*n_elem) .")
    @test size(mesh.nodes, 1) == n_nodes
    @test size(mesh.elements, 1) == n_elem
    @test dof_map.n_dofs == 2*n_edges + 4*n_faces + 6*n_elem
    @test NedelecHex.GlobalDOFs.n_edge_dofs(dof_map) == 2*n_edges
    @test NedelecHex.GlobalDOFs.n_face_dofs(dof_map) == 4*n_faces
    @test NedelecHex.GlobalDOFs.n_vol_dofs(dof_map) == 6*n_elem

    dir = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
    free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)

    # Masa definida positiva / rigidez semidefinida en los DOFs libres
    # La energía de curl es no negativa (∫|∇×u|² ≥ 0); el núcleo de K son los
    # modos gradiente espurios, por eso K solo es semidefinida.
    λM = eigvals(Symmetric(Matrix(M[free, free])))
    λK = eigvals(Symmetric(Matrix(K[free, free])))
    explica("Sobre los DOFs libres (PEC eliminados): M definida positiva, λ_min = $(minimum(λM)) > 0; K semidefinida, λ_min = $(minimum(λK)) ≳ 0.")
    @test minimum(λM) > 0
    @test minimum(λK) > -1e-8

    # Normas elementales de referencia (test de regresión)
    # Valores calculados con cuadratura de orden 4 sobre el hexaedro unitario:
    # cualquier cambio en cuadratura, W o Piola movería estas normas.
    nodes1 = [NedelecHex.ReferenceHex.NODES[k] for k in 1:8]
    Ke = NedelecHex.Assembly.stiff_matrix_element(W_REF, CURL_MONO, nodes1)
    Me = NedelecHex.Assembly.mass_matrix_element(W_REF, MONO, nodes1)
    explica("Comparación directa con la referencia (elemento unitario): ‖Kᵉ‖ = $(round(norm(Ke), digits=6)) y ‖Mᵉ‖ = $(round(norm(Me), digits=6)) frente a 556.045466 y 16.626423 .")
    @test isapprox(norm(Ke), 556.045466, rtol=1e-5)
    @test isapprox(norm(Me), 16.626423, rtol=1e-5)
end
