module Mesh
# pre-generadas en meshesHalfFilledTAP2021.mat; no hay generador con Gmsh).
# element_nodes -> +geometry/@Element/getAllCoordinates.m;
# load_hex_mesh lee los CSV que produce exportMesh.m.

using Gmsh: gmsh
using DelimitedFiles

using ..ReferenceHex: FACES

export HexMesh, PolyFace, Solid, box_mesh, load_hex_mesh, load_matlab_mesh,
    element_nodes, element_epsilon, element_inv_mu,
    edge_polyface_type, face_polyface_type, has_boundary_tags,
    tag_exterior_boundary!, tag_solids_by_plane!, edge_key, face_key

# ---------------------------------------------------------------------------
#
# y el ensamblaje solo las consulta.
#   +geometry/@PolyFace/PolyFace.m  → condición de contorno (polyFaceType)
#   +geometry/@Solid/Solid.m        → material (materialObj)
#   @Edge/@Face .polyFaceId         → 0 = interior, k = pertenece a polyFaces(k)
#   @Element   .solidId             → índice en solids
# Aquí se replica esa misma estructura.
# ---------------------------------------------------------------------------

# Condición de contorno de un conjunto de caras/aristas.
# `type` es la cadena que consulta el ensamblaje: "PerfectE" (PEC, Dirichlet),
# y deja sitio para "PMC", "Port", "ABC"… sin tocar el ensamblaje.
struct PolyFace
    id::Int
    name::String
    type::String
end

# Material de un conjunto de elementos.
struct Solid
    id::Int
    name::String
    epsilon_r::Float64
    mu_r::Float64
end

struct HexMesh
    nodes::Matrix{Float64}   # n_nodos × 3
    elements::Matrix{Int}    # n_elementos × 8  (base 1, orden de ReferenceHex.NODES)

    # --- etiquetas de material ---
    solids::Vector{Solid}
    element_solid::Vector{Int}                  # n_elem: índice en solids

    # --- etiquetas de contorno ---
    polyfaces::Vector{PolyFace}
    edge_polyface::Dict{Tuple{Int,Int},Int}     # arista (nodo_min, nodo_max) → id de polyface
    face_polyface::Dict{NTuple{4,Int},Int}      # cara (4 nodos ordenados)    → id de polyface
end

# Malla sin etiquetar: un único sólido de vacío y ninguna condición de contorno.
# Las etiquetas se ponen después (box_mesh, load_matlab_mesh,
# tag_exterior_boundary!, tag_solids_by_plane!).
function HexMesh(nodes::Matrix{Float64}, elements::Matrix{Int})
    n_elem = size(elements, 1)
    return HexMesh(nodes, elements,
        [Solid(1, "vacio", 1.0, 1.0)], ones(Int, n_elem),
        PolyFace[], Dict{Tuple{Int,Int},Int}(), Dict{NTuple{4,Int},Int}())
end

Base.show(io::IO, m::HexMesh) =
    print(io, "HexMesh($(size(m.nodes,1)) nodes, $(size(m.elements,1)) elements, " *
        "$(length(m.solids)) solids, $(length(m.polyfaces)) polyfaces)")

# --- claves canónicas de las entidades compartidas --------------------------
# Las mismas que usa build_dof_map para identificar aristas y caras entre
# elementos vecinos, así que una entidad etiquetada se reconoce desde cualquiera
# de los dos elementos que la tocan.
edge_key(g1::Int, g2::Int) = minmax(g1, g2)
face_key(gfn::NTuple{4,Int}) = NTuple{4,Int}(sort!(collect(gfn)))

# Tipo de condición de contorno de una arista/cara, "" si es interior.
_polyface_type(mesh::HexMesh, id::Int) = id == 0 ? "" : mesh.polyfaces[id].type

edge_polyface_type(mesh::HexMesh, g1::Int, g2::Int) =
    _polyface_type(mesh, get(mesh.edge_polyface, edge_key(g1, g2), 0))

face_polyface_type(mesh::HexMesh, gfn::NTuple{4,Int}) =
    _polyface_type(mesh, get(mesh.face_polyface, face_key(gfn), 0))

has_boundary_tags(mesh::HexMesh) = !isempty(mesh.polyfaces)

# --- material por elemento --------------------------------------------------
element_epsilon(mesh::HexMesh, e::Int) = mesh.solids[mesh.element_solid[e]].epsilon_r
element_epsilon(mesh::HexMesh) = [element_epsilon(mesh, e) for e in axes(mesh.elements, 1)]
element_inv_mu(mesh::HexMesh, e::Int) = 1.0 / mesh.solids[mesh.element_solid[e]].mu_r
element_inv_mu(mesh::HexMesh) = [element_inv_mu(mesh, e) for e in axes(mesh.elements, 1)]

# Los 8 nodos del elemento en coordenadas físicas, en el orden de
# ReferenceHex.NODES.
element_nodes(mesh::HexMesh, e::Int) =
    [mesh.nodes[mesh.elements[e, k], :] for k in 1:8]


# ---------------------------------------------------------------------------
# Etiquetado (preproceso). Es el único sitio del paquete donde se mira la
# geometría para decidir contorno o material: el ensamblaje solo lee etiquetas.
# ---------------------------------------------------------------------------

# Marca con una PolyFace todas las caras exteriores de la malla y las aristas
# que las bordean. El criterio es topológico, no geométrico: una cara es
# exterior si pertenece a un solo elemento. No usa coordenadas ni tolerancias,
# así que vale para cualquier forma (no solo cajas alineadas con los ejes).
function tag_exterior_boundary!(mesh::HexMesh; name="contorno", type="PerfectE")
    n_elem = size(mesh.elements, 1)

    # Cuántos elementos tocan cada cara (2 = interior, 1 = exterior).
    count = Dict{NTuple{4,Int},Int}()
    for e in 1:n_elem, face_lnodes in FACES
        k = face_key(ntuple(i -> mesh.elements[e, face_lnodes[i]], 4))
        count[k] = get(count, k, 0) + 1
    end

    id = length(mesh.polyfaces) + 1
    push!(mesh.polyfaces, PolyFace(id, name, type))

    for e in 1:n_elem, face_lnodes in FACES
        gfn = ntuple(i -> mesh.elements[e, face_lnodes[i]], 4)
        k = face_key(gfn)
        count[k] == 1 || continue
        mesh.face_polyface[k] = id
        # Las 4 aristas del borde de la cara (FACES da los nodos en orden cíclico).
        for i in 1:4
            j = i == 4 ? 1 : i + 1
            mesh.edge_polyface[edge_key(gfn[i], gfn[j])] = id
        end
    end
    return mesh
end

# Reparte los elementos en dos sólidos según a qué lado del plano `axis = cut`
# cae su centroide. El eje es un argumento explícito, no una suposición
# asumiendo z) fue exactamente eso.
# Comprueba además que ningún elemento cruce el plano, es decir, que el corte
# coincide con una capa de la malla; si no, la asignación sería ambigua.
function tag_solids_by_plane!(mesh::HexMesh; axis::Symbol=:z, cut::Float64,
    below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0),
    tol=1e-12, strict::Bool=true)

    ax = axis === :x ? 1 : axis === :y ? 2 : axis === :z ? 3 :
                                             error("axis debe ser :x, :y o :z (recibido $axis)")

    empty!(mesh.solids)
    push!(mesh.solids, Solid(1, below[1], below[2], below[3]))
    push!(mesh.solids, Solid(2, above[1], above[2], above[3]))

    n_elem = size(mesh.elements, 1)
    for e in 1:n_elem
        coord = @view mesh.nodes[mesh.elements[e, :], ax]
        lo, hi = extrema(coord)
        if lo < cut - tol && hi > cut + tol
            msg = "El elemento $e cruza el plano $axis = $cut (rango [$lo, $hi]): " *
                "el corte no cae en una capa de la malla y el material es ambiguo."
            strict ? error(msg) : @warn msg
        end
        centroid = sum(coord) / 8
        mesh.element_solid[e] = centroid < cut ? 1 : 2
    end
    return mesh
end


#genera una malla estructurada de una caja [0,Lx]×[0,Ly]×[0,Lz] con nx×ny×nz elementos hex8, usando Gmsh.
# Las 6 caras del bloque se registran como un grupo físico de Gmsh (`boundary_name`)
# y las etiquetas se leen de vuelta de la malla generada: no se deduce el contorno
# de las coordenadas.
function box_mesh(; Lx=1.0, Ly=1.0, Lz=1.0, nx=2, ny=2, nz=2,
    boundary_name="contorno", boundary_type="PerfectE")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, Lx, Ly, Lz)
    gmsh.model.occ.synchronize()

    # Asigna divisiones transfinitas a cada curva según su dirección dominante.
    lc = [Lx/nx, Ly/ny, Lz/nz]
    #bucle que recorre todas  las aristas y para cada saca la dirección y el número de puntos con Lc y divide esa arista en n_pts
    for (_, tag) in gmsh.model.getEntities(1)
        bb = gmsh.model.getBoundingBox(1, tag)
        d = [bb[4]-bb[1], bb[5]-bb[2], bb[6]-bb[3]]
        dir = argmax(d)
        n_pts = round(Int, d[dir] / lc[dir]) + 1
        gmsh.model.mesh.setTransfiniteCurve(tag, n_pts)#como cada arista tiene el mismo numero de divisiones es estructurada
    end
    for (_, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)#lo mismo para las caras
        gmsh.model.mesh.setRecombine(2, tag)
    end
    for (_, tag) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(tag)#lo mismo para el volumen
    end

    # Grupo físico con las 6 superficies del bloque = condición de contorno.
    # Es la forma estándar de etiquetar en Gmsh y el equivalente de la PolyFace
    surf_tags = [tag for (_, tag) in gmsh.model.getEntities(2)]
    pg = gmsh.model.addPhysicalGroup(2, surf_tags)
    gmsh.model.setPhysicalName(2, pg, boundary_name)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.recombine()

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    n_nodes = length(node_tags)
    nodes_mat = Matrix(reshape(coords, 3, n_nodes)')   # esto le da el formato de HexMesh

    # Gmsh usa etiquetas arbitrarias, usamos un diccionario para relacionar cada etiqueta 
    # y hacer que sean enteros consecutivas 
    tag_to_idx = Dict{Int,Int}()
    sizehint!(tag_to_idx, n_nodes)#optimiza la asignación de memoria para el diccionario
    for (i, t) in enumerate(node_tags)
        tag_to_idx[t] = i
    end

    # Elementos (hex8 = gmsh tipo 5)
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(3)
    hex_pos = findfirst(==(5), elem_types)
    isnothing(hex_pos) && error("No se encontraron elementos hex8 en la malla generada por Gmsh.")

    raw = elem_node_tags[hex_pos] #vector plano de etiquetas de nodos para todos los elementos hex8
    n_elem = length(raw) ÷ 8 #lo connvierte en una matriz de 8 x n_elem donde cada columna es un elemento y cada fila es un nodo del elemento
    conn = reshape(raw, 8, n_elem)   # etiqueta Gmsh del vértice k del elemento e = conn[k,e]
    elements = Matrix{Int}(undef, n_elem, 8)
    for e in 1:n_elem, k in 1:8
        elements[e, k] = tag_to_idx[conn[k, e]] #pasar al diccionario
    end

    # Gmsh transfinite hex8 devuelve los nodos en un orden rotado respecto al
    # convenio de ReferenceHex.NODES. Se permuta a: (0,0,0),(1,0,0),(1,1,0),(0,1,0),
    # (0,0,1),(1,0,1),(1,1,1),(0,1,1)
    for e in 1:n_elem
        elements[e, :] = elements[e, [3, 7, 8, 4, 2, 6, 5, 1]]
    end

    mesh = HexMesh(nodes_mat, elements)

    # Etiquetas de contorno: los quad4 (tipo 3) que Gmsh ha mallado sobre las
    # superficies del grupo físico. Cada quad da una cara y sus 4 aristas.
    push!(mesh.polyfaces, PolyFace(1, boundary_name, boundary_type))
    for surf in gmsh.model.getEntitiesForPhysicalGroup(2, pg)
        f_types, _, f_node_tags = gmsh.model.mesh.getElements(2, surf)
        quad_pos = findfirst(==(3), f_types)
        isnothing(quad_pos) && continue
        qraw = f_node_tags[quad_pos]
        for q in 1:(length(qraw)÷4)
            gfn = ntuple(i -> tag_to_idx[qraw[4*(q-1)+i]], 4)
            mesh.face_polyface[face_key(gfn)] = 1
            for i in 1:4
                j = i == 4 ? 1 : i + 1
                mesh.edge_polyface[edge_key(gfn[i], gfn[j])] = 1
            end
        end
    end

    gmsh.finalize()
    return mesh
end

# Reordena los 8 nodos locales de cada elemento al convenio de
# ReferenceHex.NODES a partir de su posición  dentro de la caja
# borde del elemento (asume hexaedros alineados con los ejes).
#
# referencia (ver Element.m: coordinates = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1;
# 1 0 1; 1 1 1; 0 1 1]), pero la numeración de
# meshesHalfFilledTAP2021.mat me salía diferente a veces
const _CORNER_CODE_TO_REF = Dict((0, 0, 0)=>1, (1, 0, 0)=>2, (1, 1, 0)=>3, (0, 1, 0)=>4,
    (0, 0, 1)=>5, (1, 0, 1)=>6, (1, 1, 1)=>7, (0, 1, 1)=>8)

function _canonical_node_order(nodes::Matrix{Float64}, conn::AbstractVector{Int})
    xs = @view nodes[conn, 1]
    ys = @view nodes[conn, 2]
    zs = @view nodes[conn, 3]
    xmid, ymid, zmid = (minimum(xs)+maximum(xs))/2, (minimum(ys)+maximum(ys))/2, (minimum(zs)+maximum(zs))/2
    perm = zeros(Int, 8)
    for k in 1:8
        code = (xs[k] > xmid ? 1 : 0, ys[k] > ymid ? 1 : 0, zs[k] > zmid ? 1 : 0)
        perm[_CORNER_CODE_TO_REF[code]] = k
    end
    return conn[perm]
end

function load_hex_mesh(nodes_csv::String, elems_csv::String)
    nodes = readdlm(nodes_csv, ',', Float64)
    elements = readdlm(elems_csv, ',', Int)
    for e in axes(elements, 1)
        elements[e, :] = _canonical_node_order(nodes, elements[e, :])
    end
    return HexMesh(nodes, elements)
end

#
#
# Si faltan las etiquetas de contorno se marca el exterior topológico como
# PerfectE (equivale a lo que hay en estas mallas: cavidad PEC completa) y se
# avisa. Si faltan las de material se reconstruyen los sólidos a partir de
function load_matlab_mesh(dir::String, i::Int; boundary_type="PerfectE")
    f(name) = joinpath(dir, "matlab_$(name)_$i.csv")
    mesh = load_hex_mesh(f("nodes"), f("elems"))
    n_elem = size(mesh.elements, 1)

    # Un CSV vacío cuenta como ausente (writematrix deja el fichero creado si no
    # había nada que exportar).
    hay(name) = isfile(f(name)) && filesize(f(name)) > 0

    # --- contorno ---
    if hay("polyfaces") && hay("edgetags") && hay("facetags")
        pf = readdlm(f("polyfaces"), ',')
        for r in axes(pf, 1)
            id = Int(pf[r, 1])
            id == length(mesh.polyfaces) + 1 ||
                error("$(f("polyfaces")): los id deben ser 1..n en orden (fila $r → $id).")
            push!(mesh.polyfaces, PolyFace(id, strip(string(pf[r, 3])), strip(string(pf[r, 2]))))
        end
        et = readdlm(f("edgetags"), ',', Int)
        for r in axes(et, 1)
            mesh.edge_polyface[edge_key(et[r, 1], et[r, 2])] = et[r, 3]
        end
        ft = readdlm(f("facetags"), ',', Int)
        for r in axes(ft, 1)
            mesh.face_polyface[face_key(ntuple(k -> ft[r, k], 4))] = ft[r, 5]
        end
    else
        @warn "Faltan las etiquetas de contorno de la malla $i; se marca el " *
            "exterior topológico como $boundary_type."
        tag_exterior_boundary!(mesh; name="contorno", type=boundary_type)
    end

    # --- material ---
    if isfile(f("solids")) && isfile(f("elemsolid"))
        sd = readdlm(f("solids"), ',')
        empty!(mesh.solids)
        for r in axes(sd, 1)
            push!(mesh.solids, Solid(Int(sd[r, 1]), strip(string(sd[r, 4])),
                Float64(sd[r, 2]), Float64(sd[r, 3])))
        end
        es = vec(readdlm(f("elemsolid"), ',', Int))
        length(es) == n_elem ||
            error("$(f("elemsolid")) tiene $(length(es)) valores y la malla $i tiene $n_elem elementos.")
        mesh.element_solid .= es
    elseif isfile(f("eps"))
        @warn "Faltan las etiquetas de material de la malla $i; se reconstruyen " *
            "los sólidos desde matlab_eps_$i.csv."
        eps_vec = vec(readdlm(f("eps"), ',', Float64))
        length(eps_vec) == n_elem ||
            error("$(f("eps")) tiene $(length(eps_vec)) valores y la malla $i tiene $n_elem elementos.")
        vals = sort(unique(eps_vec))
        empty!(mesh.solids)
        for (k, v) in enumerate(vals)
            push!(mesh.solids, Solid(k, "eps_r=$v", v, 1.0))
        end
        for e in 1:n_elem
            mesh.element_solid[e] = findfirst(==(eps_vec[e]), vals)
        end
    else
        error("No hay material para la malla $i: falta $(f("solids"))/$(f("elemsolid")) " *
            "y también $(f("eps")).")
    end

    return mesh
end

end # module Mesh
