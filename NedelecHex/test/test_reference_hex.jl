@testset "ReferenceHex.jl — Geometría del hexaedro de referencia" begin
    NODES        = NedelecHex.ReferenceHex.NODES
    EDGES        = NedelecHex.ReferenceHex.EDGES
    FACES        = NedelecHex.ReferenceHex.FACES
    FACE_NORMALS = NedelecHex.ReferenceHex.FACE_NORMALS

    # ── Conteos ───────────────────────────────────────────────────────────────
    @test length(NODES) == 8
    @test length(EDGES) == 12
    @test length(FACES) == 6
    @test length(FACE_NORMALS) == 6

    # ── Los 8 nodos son exactamente los vértices de [0,1]³ ──────────────────────
    @test Set(Tuple.(NODES)) == Set(Tuple(c) for c in Iterators.product(0.0:1.0, 0.0:1.0, 0.0:1.0))

    # ── Cada arista conecta exactamente 2 vértices distintos, índices válidos ──
    for (n1, n2) in EDGES
        @test n1 != n2
        @test 1 <= n1 <= 8 && 1 <= n2 <= 8
    end

    # ── Cada arista de referencia es geométricamente una arista del cubo unidad
    #    (los dos nodos difieren en exactamente una coordenada) ─────────────────
    for (n1, n2) in EDGES
        diff = NODES[n1] .- NODES[n2]
        @test count(!=(0.0), diff) == 1
    end

    # ── Cada cara tiene exactamente 4 aristas, y esas 4 aristas (como pares no
    #    ordenados) pertenecen al conjunto EDGES ────────────────────────────────
    edge_set = Set(Set.(EDGES))
    for face in FACES
        @test length(face) == 4
        for i in 1:4
            a, b = face[i], face[mod1(i + 1, 4)]
            @test Set([a, b]) in edge_set
        end
    end

    # ── Orientación local de las aristas ─────────────────────────────────────
    # La mayoría de las aristas van de nodo menor a nodo mayor, pero dos
    # (índices 3 y 7: (4,3) y (8,7)) se definen invertidas a propósito para que
    # el recorrido de las caras 1 y 2 sea un bucle consistente 1→2→3→4→1.
    # GlobalDOFs.build_dof_map detecta esta inversión por arista (comparando
    # contra los nodos globales de la malla) y permuta los DOFs correspondientes.
    forward  = [n1 < n2 for (n1, n2) in EDGES]
    reversed_idx = findall(!, forward)
    @test reversed_idx == [3, 7]
    @test count(forward) == 10

    # ── Normales de cara: vectores unitarios, uno por eje cartesiano en cada
    #    sentido ──────────────────────────────────────────────────────────────
    for n in FACE_NORMALS
        @test length(n) == 3
        @test isapprox(norm(n), 1.0)
    end
    @test Set(Tuple.(FACE_NORMALS)) == Set([
        (0.0, 0.0, -1.0), (0.0, 0.0, 1.0), (0.0, -1.0, 0.0),
        (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)
    ])
end
