@testset "Mesh.jl: Permutación de orden de nodos de Gmsh" begin
    # Mecanismo de la permutación, con un array de índices en vez de
    # coordenadas (más fácil de leer que test_jacobian.jl)
    # Mesh.jl hace: elements[e,:] = elements[e, [3,7,8,4,2,6,5,1]]
    # es decir: nueva_posición[i] = original[perm[i]].
    explica("Mecánica de la permutación: nueva_posición[i] = original[perm[i]]: con la identidad como entrada, el resultado es la propia permutación.")
    gmsh_slots = [1, 2, 3, 4, 5, 6, 7, 8]
    perm       = [3, 7, 8, 4, 2, 6, 5, 1]
    permuted   = gmsh_slots[perm]
    @test permuted == perm   # con la identidad como entrada, el resultado ES la permutación

    # ¿Es realmente necesaria? Comprobación directa contra Gmsh "en crudo"
    # Reproduce box_mesh(nx=1,ny=1,nz=1) pero SIN aplicar la permutación, para
    # verificar (no asumir) que el orden crudo de Gmsh no coincide con
    # ReferenceHex.NODES, y que aplicar la permutación lo corrige.
    using Gmsh: gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("box_test")
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    gmsh.model.occ.synchronize()
    for (_, tag) in gmsh.model.getEntities(1)
        gmsh.model.mesh.setTransfiniteCurve(tag, 2)
    end
    for (_, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    for (_, tag) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(tag)
    end
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.recombine()

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes_mat = Matrix(reshape(coords, 3, length(node_tags))')
    tag_to_idx = Dict(t => i for (i, t) in enumerate(node_tags))

    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(3)
    hex_pos = findfirst(==(5), elem_types)
    raw = elem_node_tags[hex_pos]
    conn_raw = [tag_to_idx[t] for t in raw[1:8]]
    gmsh.finalize()

    nodes_raw      = [nodes_mat[conn_raw[k], :] for k in 1:8]
    nodes_permuted = nodes_raw[perm]

    if TEST_VERBOSE
        println("\nNodos de entrada (ReferenceHex.NODES: el objetivo que se busca):")
        for k in 1:8
            println("  posición $k -> ", NedelecHex.ReferenceHex.NODES[k])
        end

        println("\nNodos de salida SIN permutar (orden crudo que devuelve Gmsh):")
        for k in 1:8
            println("  slot $k -> ", nodes_raw[k])
        end

        println("\nNodos de salida CON permutar (aplicando [3,7,8,4,2,6,5,1]):")
        for k in 1:8
            println("  posición $k -> ", nodes_permuted[k])
        end
        println()
    end

    explica("Contra Gmsh en crudo (malla transfinita real, no una suposición): el orden crudo NO coincide con ReferenceHex.NODES; tras aplicar [3,7,8,4,2,6,5,1] coincide exactamente: la permutación es necesaria y suficiente.")

    # El orden crudo de Gmsh NO coincide con ReferenceHex.NODES...
    @test nodes_raw != NedelecHex.ReferenceHex.NODES
    # ...pero tras aplicar la permutación, sí coincide exactamente.
    @test nodes_permuted == NedelecHex.ReferenceHex.NODES

    explica("Consecuencia sobre el Jacobiano: sin permutar, J es una rotación (no diagonal) que rompería el ensamblaje de DOFs; permutado, J = diag(Δx,Δy,Δz).")
    J_raw      = NedelecHex.BasisEval.compute_jacobian(nodes_raw, [0.5, 0.5, 0.5])
    J_permuted = NedelecHex.BasisEval.compute_jacobian(nodes_permuted, [0.5, 0.5, 0.5])
    @test !isapprox(J_raw, Diagonal(diag(J_raw)), atol=1e-10)         # sin permutar: NO es diagonal
    @test isapprox(J_permuted, Diagonal(diag(J_permuted)), atol=1e-10) # permutado: SÍ es diagonal

    # Robustez: la MISMA permutación fija funciona en todos los elementos
    # de una malla con varios elementos y tamaños distintos por eje
    explica("Robustez: la MISMA permutación fija vale para los 24 elementos de una malla 3×2×4 con tamaños distintos por eje (J diagonal en todos).")
    mesh   = NedelecHex.Mesh.box_mesh(nx=3, ny=2, nz=4)
    n_elem = size(mesh.elements, 1)
    @test n_elem == 24
    for e in 1:n_elem
        nodes = NedelecHex.Mesh.element_nodes(mesh, e)
        J     = NedelecHex.BasisEval.compute_jacobian(nodes, [0.5, 0.5, 0.5])
        @test isapprox(J, Diagonal(diag(J)), atol=1e-10)
    end
end
