module Mesh

using Gmsh: gmsh
using DelimitedFiles

export HexMesh, box_mesh, load_hex_mesh, element_nodes

struct HexMesh
    nodes::Matrix{Float64}   # n_nodes × 3
    elements::Matrix{Int}    # n_elements × 8  (1-based, ReferenceHex.NODES ordering)
end

Base.show(io::IO, m::HexMesh) =
    print(io, "HexMesh($(size(m.nodes,1)) nodes, $(size(m.elements,1)) elements)")

# Los 8 nodos del elemento e en coordenadas físicas, en el orden local de
# ReferenceHex.NODES.
element_nodes(mesh::HexMesh, e::Int) =
    [mesh.nodes[mesh.elements[e, k], :] for k in 1:8]

function box_mesh(; Lx=1.0, Ly=1.0, Lz=1.0, nx=2, ny=2, nz=2)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, Lx, Ly, Lz)
    gmsh.model.occ.synchronize()

    # Assign transfinite divisions to every curve by dominant direction.
    lc = [Lx/nx, Ly/ny, Lz/nz]
    for (_, tag) in gmsh.model.getEntities(1)
        bb = gmsh.model.getBoundingBox(1, tag)
        d  = [bb[4]-bb[1], bb[5]-bb[2], bb[6]-bb[3]]
        dir = argmax(d)
        n_pts = round(Int, d[dir] / lc[dir]) + 1
        gmsh.model.mesh.setTransfiniteCurve(tag, n_pts)
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

    # Nodos
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    n_nodes = length(node_tags)
    nodes_mat = Matrix(reshape(coords, 3, n_nodes)')   # n_nodes × 3

    # tag → sequential 1-based index
    tag_to_idx = Dict{Int,Int}()
    sizehint!(tag_to_idx, n_nodes)
    for (i, t) in enumerate(node_tags)
        tag_to_idx[t] = i
    end

    # Elementos (hex8 = gmsh tipo 5)
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(3)
    hex_pos = findfirst(==(5), elem_types)
    isnothing(hex_pos) && error("No hex8 elements found; check mesh settings")

    raw    = elem_node_tags[hex_pos]
    n_elem = length(raw) ÷ 8
    conn   = reshape(raw, 8, n_elem)   # 8 × n_elem
    elements = Matrix{Int}(undef, n_elem, 8)
    for e in 1:n_elem, k in 1:8
        elements[e, k] = tag_to_idx[conn[k, e]]
    end

    # Gmsh transfinite hex8 outputs nodes in a rotated order relative to
    # ReferenceHex.NODES convention. Permute to: (0,0,0),(1,0,0),(1,1,0),(0,1,0),
    # (0,0,1),(1,0,1),(1,1,1),(0,1,1) so J is diagonal on axis-aligned meshes.
    for e in 1:n_elem
        elements[e, :] = elements[e, [3, 7, 8, 4, 2, 6, 5, 1]]
    end

    gmsh.finalize()
    return HexMesh(nodes_mat, elements)
end

# Reordena los 8 nodos locales de cada elemento al convenio de
# ReferenceHex.NODES a partir de su posición geométrica dentro de la caja
# delimitadora del elemento (asume hexaedros alineados con los ejes).
#
# geometry.Element en MATLAB usa esta misma convención para el elemento de
# referencia (ver Element.m: coordinates = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1;
# 1 0 1; 1 1 1; 0 1 1]), pero la numeración local real de
# meshesHalfFilledTAP2021.mat no la respeta de forma consistente entre
# elementos (verificado: dos permutaciones distintas conviven en la misma
# malla). Sin esta corrección el Jacobiano sale rotado en vez de diagonal,
# igual que el bug de orden de Gmsh documentado en box_mesh.
const _CORNER_CODE_TO_REF = Dict((0,0,0)=>1, (1,0,0)=>2, (1,1,0)=>3, (0,1,0)=>4,
                                  (0,0,1)=>5, (1,0,1)=>6, (1,1,1)=>7, (0,1,1)=>8)

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

# Carga una HexMesh desde los CSV exportados por MATLAB (exportMesh.m) y
# canoniza el orden local de nodos de cada elemento al convenio de
# ReferenceHex.NODES (ver _canonical_node_order).
function load_hex_mesh(nodes_csv::String, elems_csv::String)
    nodes    = readdlm(nodes_csv, ',', Float64)
    elements = readdlm(elems_csv, ',', Int)
    for e in axes(elements, 1)
        elements[e, :] = _canonical_node_order(nodes, elements[e, :])
    end
    return HexMesh(nodes, elements)
end

end # module Mesh
