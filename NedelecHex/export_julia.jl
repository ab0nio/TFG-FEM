# export_julia.jl
# ---------------------------------------------------------------------------
# Genera una malla "creativa" de hexaedros estructurados con Gmsh que escribe
# la palabra "JULIA" en un estilo de píxeles (bloques).
#
# Cada letra se define como un mapa de bits 5x7; cada celda encendida es una
# caja mallada con el mismo procedimiento transfinito que usa box_mesh
# (setTransfiniteCurve/Surface/Volume + setRecombine), de modo que toda la
# malla es de hexaedros estructurados. Se exporta a un .vtk para ParaView.
#
# Uso:
#   julia --project=NedelecHex NedelecHex/export_julia.jl
#
# Salida (en NedelecHex/output/):
#   julia.vtk
# ---------------------------------------------------------------------------

using Printf
using Gmsh: gmsh

const OUT_DIR = joinpath(@__DIR__, "output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Mapa de bits de cada letra: 7 filas (fila 0 = arriba) x 5 columnas.
# '#' = celda llena, '.' = vacía.
# ---------------------------------------------------------------------------
const LETTERS = [
    # J
    ["#####",
        "....#",
        "....#",
        "....#",
        "....#",
        "#...#",
        ".###."],
    # U
    ["#...#",
        "#...#",
        "#...#",
        "#...#",
        "#...#",
        "#...#",
        ".###."],
    # L
    ["#....",
        "#....",
        "#....",
        "#....",
        "#....",
        "#....",
        "#####"],
    # I
    ["#####",
        "..#..",
        "..#..",
        "..#..",
        "..#..",
        "..#..",
        "#####"],
    # A
    [".###.",
        "#...#",
        "#...#",
        "#####",
        "#...#",
        "#...#",
        "#...#"],
]

const NROWS = 7
const NCOLS = 5
const CELL = 1.0        # lado de cada celda (cubo)
const GAP = 1.0         # separación entre letras (en celdas)
const DEPTH = 1.0       # profundidad de las letras (eje z)

# ---------------------------------------------------------------------------
# Escritura de un .vtk legacy (unstructured grid) con hexaedros (tipo 12) y un
# campo escalar por celda `letter` (1..5) para colorear en ParaView.
# ---------------------------------------------------------------------------
function write_vtk(path::String, nodes::Matrix{Float64}, elements::Matrix{Int},
    letter::Vector{Int})
    n_nodes = size(nodes, 1)
    n_elem = size(elements, 1)
    @assert size(nodes, 2) == 3
    @assert size(elements, 2) == 8
    @assert length(letter) == n_elem

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "NedelecHex JULIA (structured hexahedra)")
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
        println(io, "SCALARS letter int 1")
        println(io, "LOOKUP_TABLE default")
        for e in 1:n_elem
            println(io, letter[e])
        end
    end
    println("Escrito: ", path)
end

# ---------------------------------------------------------------------------
# Construye el modelo y genera la malla.
# ---------------------------------------------------------------------------
function build_julia()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("julia")

    # 1) Geometría: una caja por celda encendida de cada letra.
    for (li, bitmap) in enumerate(LETTERS)
        xoff = (li - 1) * (NCOLS + GAP) * CELL
        for r in 0:(NROWS-1)
            for c in 0:(NCOLS-1)
                if bitmap[r+1][c+1] == '#'
                    # fila 0 = arriba -> y crece hacia arriba.
                    gmsh.model.occ.addBox(xoff + c * CELL, (NROWS - 1 - r) * CELL,
                        0.0, CELL, CELL, DEPTH)
                end
            end
        end
    end
    gmsh.model.occ.synchronize()

    # 2) Restricciones transfinitas sobre TODO el modelo, igual que box_mesh:
    #    cada curva recibe un número de nodos según su longitud y CELL.
    for (_, tag) in gmsh.model.getEntities(1)
        bb = gmsh.model.getBoundingBox(1, tag)
        d = [bb[4] - bb[1], bb[5] - bb[2], bb[6] - bb[3]]
        dir = argmax(d)
        n_pts = max(2, round(Int, d[dir] / CELL) + 1)
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

    # 6) Asignar cada elemento a su letra por centroide.
    centroids = zeros(n_elem, 3)
    for e in 1:n_elem
        for k in 1:8
            centroids[e, :] .+= nodes_mat[elements[e, k], :]
        end
        centroids[e, :] ./= 8
    end
    letter = zeros(Int, n_elem)
    for (li, bitmap) in enumerate(LETTERS)
        xoff = (li - 1) * (NCOLS + GAP) * CELL
        for e in 1:n_elem
            cx, cy = centroids[e, 1], centroids[e, 2]
            if xoff <= cx <= xoff + NCOLS * CELL && 0.0 <= cy <= NROWS * CELL
                letter[e] = li
            end
        end
    end

    gmsh.finalize()
    return nodes_mat, elements, letter
end

nodes, elements, letter = build_julia()
write_vtk(joinpath(OUT_DIR, "julia.vtk"), nodes, elements, letter)
println("Letras: ", length(LETTERS), " | Elementos hexaédricos: ", size(elements, 1),
    " | Nodos: ", size(nodes, 1))
