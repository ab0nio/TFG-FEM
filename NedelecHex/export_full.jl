# export_full.jl
# ---------------------------------------------------------------------------
# Genera la malla de la cavidad semillena (dos materiales) con el contorno PEC
# visible, y la exporta a un .vtk para ParaView.
#
# La malla muestra TODO:
#   - los dos materiales (εᵣ = 2 dieléctrico, εᵣ = 1 vacío) → campo `material`
#   - el contorno PEC (Dirichlet) → campo `region` (1 = PEC, 0 = interior)
#
# Uso:
#   julia --project=NedelecHex NedelecHex/export_full.jl
#
# Salida (en NedelecHex/output/):
#   half_filled_<nx>x<ny>x<nz>.vtk
# ---------------------------------------------------------------------------

using DelimitedFiles
using Printf

include("src/Geometry/ReferenceHex.jl")
include("src/Geometry/Mesh.jl")
using .Mesh

const OUT_DIR = joinpath(@__DIR__, "output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Escribe un .vtk legacy con:
#   - hex8 (tipo 12): elementos de volumen, con `material` = εᵣ y `region` = 0
#   - quad (tipo 9):  caras del contorno PEC, con `material` = 0 y `region` = 1
# ---------------------------------------------------------------------------
function write_vtk_full(path::String, nodes::Matrix{Float64},
    hexes::Matrix{Int}, hex_mat::Vector{Float64},
    quads::Matrix{Int})
    n_nodes = size(nodes, 1)
    n_hex = size(hexes, 1)
    n_quad = size(quads, 1)
    n_cells = n_hex + n_quad
    @assert length(hex_mat) == n_hex

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "NedelecHex half-filled cavity with PEC boundary")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS ", n_nodes, " double")
        for i in 1:n_nodes
            @printf(io, "%.17g %.17g %.17g\n", nodes[i, 1], nodes[i, 2], nodes[i, 3])
        end

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

        println(io, "CELL_TYPES ", n_cells)
        for _ in 1:n_hex
            println(io, "12")
        end
        for _ in 1:n_quad
            println(io, "9")
        end

        # material: εᵣ de cada celda de volumen (0 para las caras PEC)
        println(io, "CELL_DATA ", n_cells)
        println(io, "SCALARS material double 1")
        println(io, "LOOKUP_TABLE default")
        for e in 1:n_hex
            println(io, hex_mat[e])
        end
        for _ in 1:n_quad
            println(io, "0")
        end

        # region: 0 interior, 1 contorno PEC
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
# ---------------------------------------------------------------------------
function extract_pec_faces(mesh::HexMesh)
    quads = Matrix{Int}(undef, 0, 4)
    for (gfn, _) in mesh.face_polyface
        quads = vcat(quads, reshape(collect(gfn), 1, 4))
    end
    return quads
end

# ---------------------------------------------------------------------------
# Genera la cavidad semillena con contorno PEC.
# ---------------------------------------------------------------------------
function export_half_filled()
    # Cavidad semillena: 0.01 × 0.001 × 0.01 m, corte en el plano medio.
    # (Lx, Ly, Lz, nx, ny, nz)
    configs = [
        (0.01, 0.001, 0.01, 2, 1, 2),   # 4  elementos
        (0.01, 0.001, 0.01, 4, 1, 4),   # 16 elementos
        (0.01, 0.001, 0.01, 8, 1, 8),   # 64 elementos
    ]
    for (Lx, Ly, Lz, nx, ny, nz) in configs
        mesh = box_mesh(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz)

        # Asigna los dos materiales por plano en Z (dieléctrico abajo, vacío
        # arriba): coincide con el modelo analítico y con la intuición.
        # below = (nombre, εᵣ, μᵣ), above = (nombre, εᵣ, μᵣ)
        tag_solids_by_plane!(mesh; axis=:z, cut=Lz / 2,
            below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0))

        # εᵣ por elemento (índice de sólido → εᵣ)
        n_elem = size(mesh.elements, 1)
        hex_mat = [mesh.solids[mesh.element_solid[e]].epsilon_r for e in 1:n_elem]

        quads = extract_pec_faces(mesh)
        base = "half_filled_$(nx)x$(ny)x$(nz)"
        write_vtk_full(joinpath(OUT_DIR, base * ".vtk"), mesh.nodes, mesh.elements,
            hex_mat, quads)
        println("  $base: ", n_elem, " elementos, ", size(quads, 1),
            " caras PEC, materiales εᵣ = ", sort(unique(hex_mat)))
    end
end

# ---------------------------------------------------------------------------
# Genera una cavidad CÚBICA (dimensiones cuadradas, no una losa plana) semillena
# con contorno PEC. Más parecida a un problema real de resonador 3D.
# ---------------------------------------------------------------------------
function export_cubic()
    # Cavidad cúbica 0.01 × 0.01 × 0.01 m, dieléctrico en la mitad inferior (Z).
    # (Lx, Ly, Lz, nx, ny, nz)
    configs = [
        (0.01, 0.01, 0.01, 4, 4, 4),     # 64   elementos
        (0.01, 0.01, 0.01, 8, 8, 8),     # 512  elementos
        (0.01, 0.01, 0.01, 12, 12, 12),  # 1728 elementos
    ]
    for (Lx, Ly, Lz, nx, ny, nz) in configs
        mesh = box_mesh(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz)

        # Dieléctrico (εᵣ=2) en la mitad inferior, vacío (εᵣ=1) arriba.
        tag_solids_by_plane!(mesh; axis=:z, cut=Lz / 2,
            below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0))

        n_elem = size(mesh.elements, 1)
        hex_mat = [mesh.solids[mesh.element_solid[e]].epsilon_r for e in 1:n_elem]

        quads = extract_pec_faces(mesh)
        base = "cubic_half_filled_$(nx)x$(ny)x$(nz)"
        write_vtk_full(joinpath(OUT_DIR, base * ".vtk"), mesh.nodes, mesh.elements,
            hex_mat, quads)
        println("  $base: ", n_elem, " elementos, ", size(quads, 1),
            " caras PEC, materiales εᵣ = ", sort(unique(hex_mat)))
    end
end

# ---------------------------------------------------------------------------
# Genera una cavidad RECTANGULAR (más larga en Z) semillena con contorno PEC.
# El corte es PARALELO al eje Z: un plano vertical x = Lx/2 que contiene al eje Z,
# dividiendo la cavidad en mitad izquierda (dieléctrico) y derecha (vacío).
# Esto coincide con la orientación de las mallas de MATLAB.
# ---------------------------------------------------------------------------
function export_rectangular()
    # Cavidad 0.01 × 0.01 × 0.02 m (doble de larga en Z), elementos cúbicos.
    # (Lx, Ly, Lz, nx, ny, nz)
    configs = [
        (0.01, 0.01, 0.02, 6, 6, 12),    # 432  elementos
        (0.01, 0.01, 0.02, 12, 12, 24),  # 3456 elementos
    ]
    for (Lx, Ly, Lz, nx, ny, nz) in configs
        mesh = box_mesh(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz)

        # Corte PARALELO al eje Z: plano vertical x = Lx/2 (contiene al eje Z).
        # below = izquierda (dieléctrico εᵣ=2), above = derecha (vacío εᵣ=1).
        tag_solids_by_plane!(mesh; axis=:x, cut=Lx / 2,
            below=("dielectrico", 2.0, 1.0), above=("vacio", 1.0, 1.0))

        n_elem = size(mesh.elements, 1)
        hex_mat = [mesh.solids[mesh.element_solid[e]].epsilon_r for e in 1:n_elem]

        quads = extract_pec_faces(mesh)
        base = "rect_half_filled_$(nx)x$(ny)x$(nz)"
        write_vtk_full(joinpath(OUT_DIR, base * ".vtk"), mesh.nodes, mesh.elements,
            hex_mat, quads)
        println("  $base: ", n_elem, " elementos, ", size(quads, 1),
            " caras PEC, materiales εᵣ = ", sort(unique(hex_mat)))
    end
end

println("=== Generando cavidad semillena con contorno PEC ===")
export_half_filled()

println("\n=== Generando cavidad CÚBICA semillena con contorno PEC ===")
export_cubic()

println("\n=== Generando cavidad RECTANGULAR semillena con contorno PEC ===")
export_rectangular()

println("\nListo. Archivos .vtk en: ", OUT_DIR)
