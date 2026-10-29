# ---------------------------------------------------------------------------
# gen_topologies.jl
#
# Genera tres mallas .vtk (formato legacy, ASCII) del MISMO cubo unidad
# [0,1]^3, una por topología:
#
#   tetraedros  -> 6 tets por cubo   (VTK cell type 10)
#   prismas     -> 2 wedges por cubo (VTK cell type 13)
#   hexaedros   -> 1 hex por cubo    (VTK cell type 12)
#
# Para visualizar en ParaView: File > Open > *.vtk > Apply.
# Para ver las aristas: Representation = Surface With Edges.
#
# Uso:
#   julia gen_topologies.jl
#
# Salida (en el mismo directorio):
#   cubo_tetra.vtk
#   cubo_prisma.vtk
#   cubo_hexa.vtk
# ---------------------------------------------------------------------------

using Printf

# ---------------------------------------------------------------------------
# Escritura de un .vtk legacy (unstructured grid).
#   nodes    : n x 3  (Float64)
#   cells    : vector de vectores de índices de nodo (base 1)
#   celltype : tipo VTK de la celda (10=tet, 13=wedge, 12=hex)
# ---------------------------------------------------------------------------
function write_vtk(path::String, nodes::Matrix{Float64},
    cells::Vector{Vector{Int}}, celltype::Int)
    n_nodes = size(nodes, 1)
    n_cells = length(cells)
    @assert size(nodes, 2) == 3

    # CELLS: por celda -> nº de nodos + índices (base 0 en VTK)
    total_ints = sum(length(c) + 1 for c in cells)

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Cubo unidad - topologia $(basename(path))")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS ", n_nodes, " double")
        for i in 1:n_nodes
            @printf(io, "%.17g %.17g %.17g\n",
                nodes[i, 1], nodes[i, 2], nodes[i, 3])
        end
        println(io, "CELLS ", n_cells, " ", total_ints)
        for c in cells
            print(io, length(c))
            for v in c
                print(io, " ", v - 1)   # VTK usa índices base 0
            end
            println(io)
        end
        println(io, "CELL_TYPES ", n_cells)
        for _ in 1:n_cells
            println(io, celltype)
        end
    end
    println("Escrito: ", path, "  (", n_cells, " celdas)")
end

# ---------------------------------------------------------------------------
# Nodos del cubo unidad (8 vértices).
# Numeración VTK para hexaedro: 0-1-2-3 cara inferior, 4-5-6-7 cara superior.
# ---------------------------------------------------------------------------
const V = [
    0.0 0.0 0.0;   # 1
    1.0 0.0 0.0;   # 2
    1.0 1.0 0.0;   # 3
    0.0 1.0 0.0;   # 4
    0.0 0.0 1.0;   # 5
    1.0 0.0 1.0;   # 6
    1.0 1.0 1.0;   # 7
    0.0 1.0 1.0;   # 8
]

# ---------------------------------------------------------------------------
# 1) HEXAEDRO: 1 celda, 8 nodos (tipo 12).
# ---------------------------------------------------------------------------
hexa = [[1, 2, 3, 4, 5, 6, 7, 8]]

# ---------------------------------------------------------------------------
# 2) PRISMAS (wedges): 2 celdas de 6 nodos (tipo 13).
#    Se corta el cubo por la diagonal de la cara inferior/superior.
#    Wedge VTK: triángulo inferior (3) + triángulo superior (3).
# ---------------------------------------------------------------------------
prismas = [
    [1, 2, 4, 5, 6, 8],   # mitad 1
    [2, 3, 4, 6, 7, 8],   # mitad 2
]

# ---------------------------------------------------------------------------
# 3) TETRAEDROS: 6 tets de 4 nodos (tipo 10).
#    Descomposición canónica del cubo en 6 tetraedros (sin solape, cubren
#    todo el volumen). Cada tet comparte la diagonal principal 1-7.
# ---------------------------------------------------------------------------
tetras = [
    [1, 2, 4, 7],
    [2, 3, 4, 7],
    [2, 3, 7, 6],
    [2, 6, 7, 5],
    [4, 7, 8, 5],
    [2, 4, 7, 5],
]

# ---------------------------------------------------------------------------
# Escritura de los tres ficheros.
# ---------------------------------------------------------------------------
write_vtk("cubo_hexa.vtk", V, hexa, 12)
write_vtk("cubo_prisma.vtk", V, prismas, 13)
write_vtk("cubo_tetra.vtk", V, tetras, 10)

println("\nListo. Abre los .vtk en ParaView (Representation = Surface With Edges).")
