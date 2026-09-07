# export_pec.jl
# ---------------------------------------------------------------------------
# Genera una malla con box_mesh y exporta un .vtk en el que se distingue el
# contorno PEC (condición de contorno de Dirichlet) del interior, para verlo
# en ParaView.
#
# Uso:
#   julia --project=NedelecHex NedelecHex/export_pec.jl
#
# Salida (en NedelecHex/output/):
#   pec_<nx>x<ny>x<nz>.vtk   malla de volumen + caras PEC etiquetadas
# ---------------------------------------------------------------------------

using DelimitedFiles
using Printf

include("src/Geometry/ReferenceHex.jl")
include("src/Geometry/Mesh.jl")
using .Mesh
using .ReferenceHex: FACES

const OUT_DIR = joinpath(@__DIR__, "output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Escribe un .vtk legacy con dos tipos de celda:
#   - hex8 (tipo 12): elementos de volumen, region = 0
#   - quad (tipo 9):  caras del contorno PEC, region = 1
# El campo escalar `region` permite colorear el contorno en ParaView.
# ---------------------------------------------------------------------------
function write_vtk_pec(path::String, nodes::Matrix{Float64},
    hexes::Matrix{Int}, quads::Matrix{Int})
    n_nodes = size(nodes, 1)
    n_hex = size(hexes, 1)
    n_quad = size(quads, 1)
    n_cells = n_hex + n_quad

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "NedelecHex mesh with PEC boundary")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS ", n_nodes, " double")
        for i in 1:n_nodes
            @printf(io, "%.17g %.17g %.17g\n", nodes[i, 1], nodes[i, 2], nodes[i, 3])
        end

        # CELLS: cada hex = 9 enteros, cada quad = 5 enteros
        println(io, "CELLS ", n_cells, " ", n_hex * 9 + n_quad * 5)
        for e in 1:n_hex
            @printf(io, "8")
            for k in 1:8
                @printf(io, " %d", hexes[e, k] - 1)
            end
            println(io)
        end
        for q in 1:n_quad
            @printf(io, "4")
            for k in 1:4
                @printf(io, " %d", quads[q, k] - 1)
            end
            println(io)
        end

        # CELL_TYPES: 12 = hexaedro, 9 = quad
        println(io, "CELL_TYPES ", n_cells)
        for _ in 1:n_hex
            println(io, "12")
        end
        for _ in 1:n_quad
            println(io, "9")
        end

        # region: 0 interior, 1 contorno PEC
        println(io, "CELL_DATA ", n_cells)
        println(io, "SCALARS region int 1")
        println(io, "LOOKUP_TABLE default")
        for _ in 1:n_hex
            println(io, "0")
        end
        for _ in 1:n_quad
            println(io, "1")
        end
    end
    println("Escrito: ", path)
end

# ---------------------------------------------------------------------------
# Extrae las caras del contorno PEC de una HexMesh.
# box_mesh etiqueta cada cara de contorno en mesh.face_polyface (clave = 4 nodos).
# ---------------------------------------------------------------------------
function extract_pec_faces(mesh::HexMesh)
    quads = Matrix{Int}(undef, 0, 4)
    for (gfn, _) in mesh.face_polyface
        # gfn es una tupla de 4 nodos (ordenados por face_key, es decir, orden
        # canónico). Para el quad de VTK da igual el orden de los 4 vértices.
        quads = vcat(quads, reshape(collect(gfn), 1, 4))
    end
    return quads
end

# ---------------------------------------------------------------------------
# Genera la malla y exporta el .vtk con el contorno PEC.
# ---------------------------------------------------------------------------
function export_pec()
    # Cavidad semillena típica: 0.01 × 0.001 × 0.01 m.
    # (Lx, Ly, Lz, nx, ny, nz)
    configs = [
        (0.01, 0.001, 0.01, 2, 1, 2),   # 4  elementos
        (0.01, 0.001, 0.01, 4, 1, 4),   # 16 elementos
        (0.01, 0.001, 0.01, 8, 1, 8),   # 64 elementos
    ]
    for (Lx, Ly, Lz, nx, ny, nz) in configs
        mesh = box_mesh(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz)
        quads = extract_pec_faces(mesh)
        base = "pec_$(nx)x$(ny)x$(nz)"
        write_vtk_pec(joinpath(OUT_DIR, base * ".vtk"), mesh.nodes, mesh.elements, quads)
        println("  contorno PEC: ", size(quads, 1), " caras cuadrangulares")
    end
end

println("=== Generando mallas con contorno PEC visible ===")
export_pec()
println("\nListo. Archivos .vtk en: ", OUT_DIR)
