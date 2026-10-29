#  Prueba Assembly.jl y GlobalDOFs.jl.
#   - Ensamblaje global,1 → 2 → 8 elementos (incluye la dimensión del núcleo de
#     K en cada tamaño: 26 / 44 / 124 modos gradiente).
#   - comparación de K/M globales MATLAB vs Julia .
using DelimitedFiles

# Malla estructurada + mapa de DOFs + matrices globales en un paso
# (W_REF, MONO, CURL_MONO son globales de runtests.jl).
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
    # (8 vértices + 12 puntos medios de arista + 6 centros de cara + 1 centro de volumen = 27 nodos), 
    #menos la función constante (gradiente nulo) → 27 − 1 = 26 modos.
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
    @test isempty(intersect(Set(conn1[49:54]), Set(conn2[49:54])))        # volumen:0

    # Acumulación COO en DOFs compartidos
    # El ensamblaje suma las contribuciones de ambos elementos cuando comparten
    # índice global: K_global[g,h] = K1[loc1] + K2[loc2].
    nidx1, nidx2 = mesh.elements[1, :], mesh.elements[2, :]
    nodes1, nodes2 = NedelecHex.Mesh.element_nodes(mesh, 1), NedelecHex.Mesh.element_nodes(mesh, 2)
    K1 = NedelecHex.Assembly.stiff_matrix_element(W_REF, CURL_MONO, nodes1)
    K2 = NedelecHex.Assembly.stiff_matrix_element(W_REF, CURL_MONO, nodes2)

    explica("para cada pareja de DOFs compartidos (g,h), K_global debe ser (k1 + k2).")
    for g in shared, h in shared
        loc1g, loc1h = findfirst(==(g), conn1), findfirst(==(h), conn1)
        loc2g, loc2h = findfirst(==(g), conn2), findfirst(==(h), conn2)
        @test isapprox(K[g, h], K1[loc1g, loc1h] + K2[loc2g, loc2h], atol=1e-10)
    end

    # Orientación consistente de las 4 aristas de la cara compartida
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

    # Dimensión del núcleo de K (modos gradiente espurios), ahora en 2 elementos
    # La base nodal Q2 conjunta tiene (2·2+1)(2·1+1)(2·1+1) = 5·3·3 = 45 nodos;
    # menos la función constante → 45 − 1 = 44 modos con curl ≡ 0.
    λK = eigvals(Symmetric(Matrix(K)))
    n_null = count(x -> abs(x) < 1e-6, λK)
    explica("Dimensión del núcleo de K (2 elem): $(n_null) autovalores ≈ 0; nodos Q2 = 5·3·3 = 45, menos la constante → 44.")
    @test n_null == 44
end

@testset "Assembly: malla 2×2×2" begin
    mesh, dof_map, K, M = assemble_box(2, 2, 2)

    # Conteo  de DOFs para una malla estructurada nx=ny=nz=2
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

    # Dimensión del núcleo de K (modos gradiente espurios), malla completa 2×2×2
    # La base nodal Q2 tiene (2·2+1)³ = 5³ = 125 nodos; menos la función
    # constante → 125 − 1 = 124 modos con curl ≡ 0.
    λK_full = eigvals(Symmetric(Matrix(K)))
    n_null = count(x -> abs(x) < 1e-6, λK_full)
    explica("Dimensión del núcleo de K (8 elem): $(n_null) autovalores ≈ 0; nodos Q2 = 5³ = 125, menos la constante → 124.")
    @test n_null == 124

    dir = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
    free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)

    # Masa definida positiva / rigidez semidefinida en los DOFs libres.
    λM = eigvals(Symmetric(Matrix(M[free, free])))
    λK = eigvals(Symmetric(Matrix(K[free, free])))
    explica("Sobre los DOFs libres (PEC eliminados): M definida positiva, λ_min = $(minimum(λM)) > 0; K semidefinida, λ_min = $(minimum(λK)) ≳ 0.")
    @test minimum(λM) > 0
    @test minimum(λK) > -1e-8
end

@testset "Assembly: dimensión del núcleo de K en mallas adicionales" begin
    # Verifica la fórmula cerrada del núcleo de K (modos gradiente espurios):
    # nodos Q2 = (2·nx+1)(2·ny+1)(2·nz+1), menos la función constante → nodos − 1.
    # Mallas nuevas añadidas a la tabla de la memoria (tab:nucleo).
    casos = [
        (2, 2, 1, 74),   # (2·2+1)(2·2+1)(2·1+1) = 5·5·3 = 75 → 74
        (3, 3, 3, 342),  # 7³ = 343 → 342
        (4, 4, 4, 728),  # 9³ = 729 → 728
    ]
    for (nx, ny, nz, n_null_expected) in casos
        mesh, dof_map, K, _ = assemble_box(nx, ny, nz)
        λK = eigvals(Symmetric(Matrix(K)))
        n_null = count(x -> abs(x) < 1e-6, λK)
        explica("Dimensión del núcleo de K (malla $(nx)×$(ny)×$(nz)): $(n_null) autovalores ≈ 0; esperado $(n_null_expected).")
        @test n_null == n_null_expected
    end
end

# Comparación de las matrices globales K y M ensambladas: MATLAB vs Julia.
# Requiere los CSV que produce matlab/exportGlobalMatrices.m (matlab_K_N.csv,
# matlab_M_N.csv, matlab_free_N.csv).
@testset "Assembly: comparación K/M globales MATLAB vs Julia (" begin
    W, mono, curl_mono = W_REF, MONO, CURL_MONO
    mesh_dir = joinpath(@__DIR__, "data", "meshes")

    #Función auxiliar para cargar una matriz dispersa en formato COO exportada por MATLAB
    load_coo(path) = (A=readdlm(path, ','); sparse(Int.(A[:, 1]), Int.(A[:, 2]), Float64.(A[:, 3])))
    # Función auxiliar para calcular la diferencia relativa entre dos valores 
    reldiff(x, y) = abs(x - y) / max(abs(x), abs(y), eps())

    # Los 4 ficheros que hace falta tener de MATLAB para comparar la malla i.
    have(i) = all(isfile(joinpath(mesh_dir, "matlab_$(f)_$i.csv")) for f in ("K", "M", "free", "eps"))

    any_data = any(have(i) for i in 1:5)
    if !any_data
        explica(" faltan output/matlab_{K,M,free,eps}_N.csv. se generan con matlab/exportGlobalMatrices.m y matlab/exportElementEpsilon.m.")
        @test_skip any_data
    else
        for i in 1:5
            have(i) || continue #Si falta algún fichero, pasa a la siguiente malla

            K_ml = load_coo(joinpath(mesh_dir, "matlab_K_$i.csv"))
            M_ml = load_coo(joinpath(mesh_dir, "matlab_M_$i.csv"))
            free_ml = Int.(vec(readdlm(joinpath(mesh_dir, "matlab_free_$i.csv"), ',')))#índices de DOFs libres en MATLAB 
            Kf_ml, Mf_ml = K_ml[free_ml, free_ml], M_ml[free_ml, free_ml]

            #carga la malla de MATLAB con sus etiquetas (sólidos y contorno) para ensamblar en Julia
            mesh = NedelecHex.Mesh.load_matlab_mesh(mesh_dir, i)
            dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)

            # εᵣ y PEC salen de las etiquetas de MATLAB, no de las coordenadas.
            K_jl, M_jl = NedelecHex.Assembly.assemble_global_matrices(
                mesh, dof_map, W, mono, curl_mono)

            free_jl = NedelecHex.Assembly.free_dofs(dof_map.n_dofs,
                NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map))

            Kf_jl, Mf_jl = K_jl[free_jl, free_jl], M_jl[free_jl, free_jl]

            @printf("\n-- malla %d --\n", i)
            @printf("  n_free   MATLAB=%d  Julia=%d\n", length(free_ml), length(free_jl))
            @printf("  ||Kf||_F MATLAB=%.6e  Julia=%.6e  difrel=%.2e\n",
                norm(Kf_ml), norm(Kf_jl), reldiff(norm(Kf_ml), norm(Kf_jl)))
            @printf("  ||Mf||_F MATLAB=%.6e  Julia=%.6e  difrel=%.2e\n",
                norm(Mf_ml), norm(Mf_jl), reldiff(norm(Mf_ml), norm(Mf_jl)))
            @printf("  tr(Kf)   MATLAB=%.6e  Julia=%.6e  difrel=%.2e\n",
                tr(Kf_ml), tr(Kf_jl), reldiff(tr(Kf_ml), tr(Kf_jl)))
            @printf("  tr(Mf)   MATLAB=%.6e  Julia=%.6e  difrel=%.2e\n",
                tr(Mf_ml), tr(Mf_jl), reldiff(tr(Mf_ml), tr(Mf_jl)))

            # Mismo problema discreto → mismas invariantes. 
            @test length(free_ml) == length(free_jl)#mismas incógnitas tras imponer PEC

            #Kf y Mf son matrices K y M restringidos a los grados de libertad libres (tras eliminar los DOFs de contorno).
            #COmprobamos que la norma y la traza sea igual entre MATLAB y Julia 
            @test reldiff(norm(Kf_ml), norm(Kf_jl)) < 1e-10
            @test reldiff(norm(Mf_ml), norm(Mf_jl)) < 1e-10
            @test reldiff(tr(Kf_ml), tr(Kf_jl)) < 1e-10
            @test reldiff(tr(Mf_ml), tr(Mf_jl)) < 1e-10
        end
    end
end
