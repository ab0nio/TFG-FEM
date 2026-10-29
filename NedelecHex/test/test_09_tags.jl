# Prueba el etiquetado de la malla: PolyFace (condición de contorno) y Solid
# (material), y que el ensamblaje los consulta en vez de deducirlos de las
# coordenadas.
#
# MATLAB: +geometry/@PolyFace/PolyFace.m, +geometry/@Solid/Solid.m y
# +assembly/imposeLocalDirichletMatrix.m (checkDirichlet compara polyFaceType
# con "PerfectE"; no mira ninguna coordenada).
using DelimitedFiles

using NedelecHex.Mesh: HexMesh, PolyFace, box_mesh, load_matlab_mesh,
    tag_exterior_boundary!, tag_solids_by_plane!, element_epsilon,
    has_boundary_tags, edge_key, face_key
using NedelecHex.GlobalDOFs: build_dof_map, dirichlet_dofs
using NedelecHex.ReferenceHex: EDGES, FACES

# Criterio geométrico que se usaba antes de las etiquetas: una arista/cara es
# Dirichlet si todos sus nodos comparten un plano frontera de la caja. Se
# conserva aquí SOLO como referencia para comprobar que el nuevo mecanismo da
# exactamente el mismo resultado en el caso donde el viejo era válido.
function dirichlet_dofs_geometrico(mesh, dof_map; tol=1e-10)
    coords = mesh.nodes
    xmin, xmax = extrema(coords[:, 1])
    ymin, ymax = extrema(coords[:, 2])
    zmin, zmax = extrema(coords[:, 3])
    on_plane = falses(size(coords, 1), 6)
    for n in axes(coords, 1)
        on_plane[n, 1] = abs(coords[n, 1] - xmin) < tol
        on_plane[n, 2] = abs(coords[n, 1] - xmax) < tol
        on_plane[n, 3] = abs(coords[n, 2] - ymin) < tol
        on_plane[n, 4] = abs(coords[n, 2] - ymax) < tol
        on_plane[n, 5] = abs(coords[n, 3] - zmin) < tol
        on_plane[n, 6] = abs(coords[n, 3] - zmax) < tol
    end

    dir = Set{Int}()
    for e in axes(mesh.elements, 1)
        gn = mesh.elements[e, :]
        dofs = dof_map.connectivity[e, :]
        for (k, (ln1, ln2)) in enumerate(EDGES)
            if any(on_plane[gn[ln1], :] .& on_plane[gn[ln2], :])
                push!(dir, dofs[2*(k-1)+1], dofs[2*(k-1)+2])
            end
        end
        for (f, face_lnodes) in enumerate(FACES)
            if any(reduce(.&, on_plane[gn[ln], :] for ln in face_lnodes))
                for i in 1:4
                    push!(dir, dofs[24+4*(f-1)+i])
                end
            end
        end
    end
    return sort!(collect(dir))
end

@testset "Etiquetas: box_mesh trae el contorno puesto por Gmsh" begin
    nx, ny, nz = 3, 2, 4
    mesh = box_mesh(nx=nx, ny=ny, nz=nz)

    @test has_boundary_tags(mesh)
    @test length(mesh.polyfaces) == 1
    @test mesh.polyfaces[1].type == "PerfectE"

    # Caras exteriores de una malla estructurada nx×ny×nz: las 6 paredes.
    n_caras = 2 * (nx*ny + ny*nz + nx*nz)
    explica("box_mesh registra las 6 superficies del bloque como grupo físico de Gmsh
    y lee de vuelta los quad4 mallados sobre ellas: $(length(mesh.face_polyface)) caras
    etiquetadas, esperadas 2·(nx·ny + ny·nz + nx·nz) = $n_caras.")
    @test length(mesh.face_polyface) == n_caras

    # Aristas de la superficie: por Euler sobre la superficie cerrada
    # (V − A + C = 2) → A = V + C − 2, con V = nodos totales menos interiores.
    n_nodos_sup = (nx+1)*(ny+1)*(nz+1) - (nx-1)*(ny-1)*(nz-1)
    explica("Aristas etiquetadas: $(length(mesh.edge_polyface)); por Euler sobre la
    superficie cerrada, V − A + C = 2 con V = $n_nodos_sup nodos de superficie y
    C = $n_caras caras → A = $(n_nodos_sup + n_caras - 2).")
    @test length(mesh.edge_polyface) == n_caras + n_nodos_sup - 2
end

@testset "Etiquetas: el Dirichlet por etiqueta = el geométrico (no regresión)" begin
    # Donde el criterio viejo era correcto (caja PEC completa alineada con los
    # ejes) el nuevo tiene que dar EXACTAMENTE el mismo conjunto de DOFs.
    for (nx, ny, nz) in ((1, 1, 1), (2, 1, 1), (2, 2, 2), (3, 2, 4))
        mesh = box_mesh(nx=nx, ny=ny, nz=nz)
        dof_map = build_dof_map(mesh)
        dir_tag = dirichlet_dofs(mesh, dof_map)
        dir_geo = dirichlet_dofs_geometrico(mesh, dof_map)
        @test dir_tag == dir_geo
    end
    explica("Malla $(("1×1×1, 2×1×1, 2×2×2, 3×2×4")): el conjunto de DOFs Dirichlet
    por etiqueta coincide entero con el que daba el criterio de planos. El cambio es de
    procedencia del dato (etiqueta vs coordenadas), no de resultado.")

    # tag_exterior_boundary! es topológico (una cara es exterior si la toca un
    # solo elemento): sin coordenadas ni tolerancia, y da lo mismo.
    mesh = box_mesh(nx=3, ny=2, nz=4)
    sin_tags = HexMesh(mesh.nodes, mesh.elements)
    tag_exterior_boundary!(sin_tags)
    dof_map = build_dof_map(sin_tags)
    @test dirichlet_dofs(sin_tags, dof_map) == dirichlet_dofs_geometrico(sin_tags, dof_map)
end

@testset "Etiquetas: sin etiquetas no se adivina, se falla" begin
    mesh = box_mesh(nx=2, ny=2, nz=2)
    sin_tags = HexMesh(mesh.nodes, mesh.elements)   # constructor sin etiquetar
    dof_map = build_dof_map(sin_tags)

    explica("Una malla sin etiquetas de contorno no tiene respuesta a '¿qué es PEC?'.
    Antes se contestaba mirando la caja envolvente, que acierta solo si el problema es
    una caja PEC entera; ahora se lanza un error en vez de suponerlo en silencio.")
    @test !has_boundary_tags(sin_tags)
    @test_throws ErrorException dirichlet_dofs(sin_tags, dof_map)
end

@testset "Etiquetas: PEC parcial (lo que el criterio geométrico no sabe hacer)" begin
    # Solo la pared x = 0 es metálica; las otras cinco, no. Es el caso que
    # motivó el cambio: la geometría de la caja es idéntica, cambia la física.
    mesh = box_mesh(nx=2, ny=2, nz=2)
    empty!(mesh.polyfaces)
    empty!(mesh.face_polyface)
    empty!(mesh.edge_polyface)
    push!(mesh.polyfaces, PolyFace(1, "pared_x0", "PerfectE"))

    xmin = minimum(@view mesh.nodes[:, 1])
    en_x0(n) = abs(mesh.nodes[n, 1] - xmin) < 1e-10
    for e in axes(mesh.elements, 1), face_lnodes in FACES
        gfn = ntuple(i -> mesh.elements[e, face_lnodes[i]], 4)
        all(en_x0, gfn) || continue
        mesh.face_polyface[face_key(gfn)] = 1
        for i in 1:4
            j = i == 4 ? 1 : i + 1
            mesh.edge_polyface[edge_key(gfn[i], gfn[j])] = 1
        end
    end

    dof_map = build_dof_map(mesh)
    dir = dirichlet_dofs(mesh, dof_map)

    # La pared x=0 de una malla 2×2×2: 2×2 = 4 caras (4 DOFs cada una) y
    # 2·(3·2) = 12 aristas (2 DOFs cada una) → 4·4 + 12·2 = 40 DOFs.
    explica("PEC solo en x=0: 4 caras × 4 DOFs + 12 aristas × 2 DOFs = 40 DOFs
    Dirichlet, frente a los $(length(dirichlet_dofs_geometrico(mesh, dof_map))) que
    impondría el criterio de la caja envolvente (las 6 paredes). Misma malla, física
    distinta: la diferencia está en la etiqueta.")
    @test length(dir) == 40
    @test length(dir) < length(dirichlet_dofs_geometrico(mesh, dof_map))
end

@testset "Etiquetas: material por sólido y el eje del corte" begin
    # El εᵣ ya no se decide dentro del ensamblaje: se etiqueta antes, con el eje
    # del corte explícito.
    mesh_z = box_mesh(nx=2, ny=2, nz=2)
    tag_solids_by_plane!(mesh_z; axis=:z, cut=0.5,
        below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0))
    mesh_x = box_mesh(nx=2, ny=2, nz=2)
    tag_solids_by_plane!(mesh_x; axis=:x, cut=0.5,
        below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0))

    eps_z, eps_x = element_epsilon(mesh_z), element_epsilon(mesh_x)
    @test count(==(2.0), eps_z) == 4
    @test count(==(2.0), eps_x) == 4

    # Dos cortes simétricos de una caja coinciden en la mitad de los elementos:
    # es exactamente el 50% que medí entre el heurístico en z y el εᵣ real de
    # las mallas de MATLAB (cortadas por x). Un 50% de acierto es el azar, no
    # una aproximación buena.
    coincidencias = count(eps_z .== eps_x)
    explica("Mismo dominio, mismo corte a la mitad, distinto eje: las dos asignaciones
    coinciden en $coincidencias de 8 elementos (50%). Por eso el heurístico por centroide
    en z 'casi' funcionaba sobre mallas cortadas por x sin dar ningún error.")
    @test coincidencias == 4

    # Un corte que no cae en una capa de la malla deja elementos a caballo: el
    # material sería ambiguo y se avisa en vez de repartirlo por el centroide.
    mesh_mal = box_mesh(nx=2, ny=2, nz=1)
    @test_throws ErrorException tag_solids_by_plane!(mesh_mal; axis=:z, cut=0.5)
end

@testset "Etiquetas: las mallas de MATLAB traen su propio material" begin
    mesh_dir = joinpath(@__DIR__, "data", "meshes")
    have(i) = isfile(joinpath(mesh_dir, "matlab_nodes_$i.csv")) &&
        isfile(joinpath(mesh_dir, "matlab_eps_$i.csv"))

    if !any(have(i) for i in 1:5)
        explica("faltan test/data/meshes/matlab_{nodes,eps}_N.csv; se generan con " *
            "scripts/generate_material_tags.jl.")
        @test_skip any(have(i) for i in 1:5)
    else
        for i in 1:5
            have(i) || continue
            mesh = load_matlab_mesh(mesh_dir, i)
            eps_jl = element_epsilon(mesh)
            eps_ml = vec(readdlm(joinpath(mesh_dir, "matlab_eps_$i.csv"), ',', Float64))
            @test eps_jl == eps_ml

            # Con las polyFaces reales de la malla, el conjunto de DOFs
            # Dirichlet es el MISMO que daba el criterio de la caja envolvente:
            # en esta cavidad (PEC en las 6 paredes) el criterio viejo acertaba,
            # solo que por casualidad del problema.
            pf = joinpath(mesh_dir, "matlab_polyfaces_$i.csv")
            if isfile(pf) && filesize(pf) > 0
                dof_map = build_dof_map(mesh)
                @test dirichlet_dofs(mesh, dof_map) == dirichlet_dofs_geometrico(mesh, dof_map)
            end

            # Comprobación del eje real de la losa en estas mallas.
            n_elem = size(mesh.elements, 1)
            centro(e, ax) = sum(@view mesh.nodes[mesh.elements[e, :], ax]) / 8
            acierto_x = count(e -> (centro(e, 1) > 0) == (eps_ml[e] == 2.0), 1:n_elem)
            acierto_z = count(e -> (centro(e, 3) < 0) == (eps_ml[e] == 2.0), 1:n_elem)
            @printf("  malla %d: %d elem, criterio en x acierta %d, en z %d\n",
                i, n_elem, acierto_x, acierto_z)
            @test acierto_x == n_elem          # la losa está cortada por x
            @test acierto_z == n_elem ÷ 2      # el criterio en z acierta la mitad (azar)
        end
        explica("El εᵣ que Julia usa ahora es el de la etiqueta, idéntico al de MATLAB
        elemento a elemento. El criterio en x acierta el 100% y el de z el 50%: la losa
        de estas mallas está cortada perpendicular a x, no a z.")
    end
end
