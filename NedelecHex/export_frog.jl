# export_frog.jl
# ---------------------------------------------------------------------------
# Genera una malla "creativa" de hexaedros estructurados con Gmsh: una rana en
# cuclillas mirando de frente (cuerpo ancho, cabeza arriba con dos ojos,
# patas delanteras abajo).
#
# Cada parte (cuerpo, cabeza, ojos, patas, pies, boca) es una caja mallada con
# el mismo procedimiento transfinito que usa box_mesh
# (setTransfiniteCurve/Surface/Volume + setRecombine), de modo que toda la
# malla es de hexaedros estructurados. Se exporta a un .vtk para ParaView.
#
# Uso:
#   julia --project=NedelecHex NedelecHex/export_frog.jl
#
# Salida (en NedelecHex/output/):
#   frog.vtk
# ---------------------------------------------------------------------------

using Printf
using Gmsh: gmsh

const OUT_DIR = joinpath(@__DIR__, "output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Partes de la rana: (x0, x1, y0, y1, z0, z1, parte)
#   x = lateral, y = arriba, z = profundidad (la rana mira a +z, al espectador).
#   parte: 1 cuerpo, 2 cabeza, 3 ojos, 4 patas, 5 pies, 6 boca
# ---------------------------------------------------------------------------
const PARTS = [
    # Cuerpo (ancho y en cuclillas)
    (-2.0, 2.0, 0.0, 1.5, -1.5, 0.0, 1),
    # Cabeza (ancha, sobre el cuerpo)
    (-1.6, 1.6, 1.5, 2.4, -1.1, 0.0, 2),
    # Ojos (dos, salientes encima de la cabeza)
    (-1.3, -0.35, 2.4, 3.0, -0.7, 0.0, 3),
    (0.35, 1.3, 2.4, 3.0, -0.7, 0.0, 3),
    # Patas delanteras (a los lados, abajo)
    (-2.5, -2.0, 0.0, 1.0, -1.0, 0.0, 4),
    (2.0, 2.5, 0.0, 1.0, -1.0, 0.0, 4),
    # Pies (planos y anchos, delante abajo)
    (-2.5, -1.0, 0.0, 0.35, -0.3, 0.5, 5),
    (1.0, 2.5, 0.0, 0.35, -0.3, 0.5, 5),
    # Boca / hocico (saliente en el frente de la cabeza)
    (-1.4, 1.4, 1.7, 2.0, 0.0, 0.3, 6),
]

# Tamaño de elemento (lado de cada celda hexaédrica). Cuanto menor, más fina.
const LC = 0.15

# ---------------------------------------------------------------------------
# Escritura de un .vtk legacy (unstructured grid) con hexaedros (tipo 12) y un
# campo escalar por celda `part` (1..6) para colorear en ParaView.
# ---------------------------------------------------------------------------
function write_vtk(path::String, nodes::Matrix{Float64}, elements::Matrix{Int},
    part::Vector{Int})
    n_nodes = size(nodes, 1)
    n_elem = size(elements, 1)
    @assert size(nodes, 2) == 3
    @assert size(elements, 2) == 8
    @assert length(part) == n_elem

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "NedelecHex frog (structured hexahedra)")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS ", n_nodes, " double")
        for i in 1:n_nodes
            @printf(io, "%.17g %.17g %.17g\n", nodes[i, 1], nodes[i, 2], nodes[i, 3])
        end
        println(io, "CELLS ", n_elem, " ", n_elem * 9)
        for e in 1:n_elem
            @printf(io, "8")
            for k in 1:8
                @printf(io, " %d", elements[e, k] - 1)   # VTK usa índices base 0
            end
            println(io)
        end
        println(io, "CELL_TYPES ", n_elem)
        for _ in 1:n_elem
            println(io, "12")   # hexaedro
        end
        println(io, "CELL_DATA ", n_elem)
        println(io, "SCALARS part int 1")
        println(io, "LOOKUP_TABLE default")
        for e in 1:n_elem
            println(io, part[e])
        end
    end
    println("Escrito: ", path)
end

# ---------------------------------------------------------------------------
# Construye el modelo y genera la malla.
# ---------------------------------------------------------------------------
function build_frog()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("frog")

    # 1) Geometría: una caja por parte.
    for (x0, x1, y0, y1, z0, z1, _) in PARTS
        gmsh.model.occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)
    end
    gmsh.model.occ.synchronize()

    # 2) Restricciones transfinitas sobre TODO el modelo, igual que box_mesh:
    #    cada curva recibe un número de nodos según su longitud y LC.
    for (_, tag) in gmsh.model.getEntities(1)
        bb = gmsh.model.getBoundingBox(1, tag)
        d = [bb[4] - bb[1], bb[5] - bb[2], bb[6] - bb[3]]
        dir = argmax(d)
        n_pts = max(2, round(Int, d[dir] / LC) + 1)
        gmsh.model.mesh.setTransfiniteCurve(tag, n_pts)
    end
    for (_, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    for (_, tag) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(tag)
    end

    # 3) Generar y recombinar para obtener hexaedros.
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.recombine()

    # 4) Leer nodos.
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    n_nodes = length(node_tags)
    nodes_mat = Matrix(reshape(coords, 3, n_nodes)')
    tag_to_idx = Dict{Int,Int}()
    for (i, t) in enumerate(node_tags)
        tag_to_idx[t] = i
    end

    # 5) Leer hexaedros (tipo 5).
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(3)
    hex_pos = findfirst(==(5), elem_types)
    isnothing(hex_pos) && error("No se encontraron elementos hex8.")
    raw = elem_node_tags[hex_pos]
    n_elem = length(raw) ÷ 8
    conn = reshape(raw, 8, n_elem)
    elements = Matrix{Int}(undef, n_elem, 8)
    for e in 1:n_elem, k in 1:8
        elements[e, k] = tag_to_idx[conn[k, e]]
    end

    # 6) Asignar cada elemento a su parte por centroide.
    centroids = zeros(n_elem, 3)
    for e in 1:n_elem
        for k in 1:8
            centroids[e, :] .+= nodes_mat[elements[e, k], :]
        end
        centroids[e, :] ./= 8
    end
    part = zeros(Int, n_elem)
    for (x0, x1, y0, y1, z0, z1, p) in PARTS
        for e in 1:n_elem
            cx, cy, cz = centroids[e, 1], centroids[e, 2], centroids[e, 3]
            if x0 <= cx <= x1 && y0 <= cy <= y1 && z0 <= cz <= z1
                part[e] = p
            end
        end
    end

    gmsh.finalize()
    return nodes_mat, elements, part
end

nodes, elements, part = build_frog()
write_vtk(joinpath(OUT_DIR, "frog.vtk"), nodes, elements, part)
println("Partes: ", length(PARTS), " | Elementos hexaédricos: ", size(elements, 1),
    " | Nodos: ", size(nodes, 1))
