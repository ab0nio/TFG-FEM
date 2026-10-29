# Convierte las mallas de referencia (CSV) a formato Gmsh .msh 2.2 y exporta las
# mallas generadas con box_mesh, guardándolo todo en output/.
#
# Uso (desde NedelecHex/):
#   julia --project=. scripts/export_msh.jl
#
# Salida (en output/):
#   matlab_mesh_1.msh ... matlab_mesh_5.msh   (mallas de referencia convertidas)
#   box_mesh_<nx>x<ny>x<nz>.msh               (mallas generadas con box_mesh)

using DelimitedFiles
using Printf

# Incluye el módulo Mesh para poder usar box_mesh y HexMesh.
include("src/Geometry/ReferenceHex.jl")
include("src/Geometry/Mesh.jl")
using .Mesh

const OUT_DIR = joinpath(@__DIR__, "output")
isdir(OUT_DIR) || mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Escritura de un .msh 2.2 a partir de nodos, elementos y etiqueta física.
#
# Formato .msh 2.2 (hex8 = tipo 5):
#   $MeshFormat
#   2.2 0 8
#   $EndMeshFormat
#   $Nodes
#   <n_nodos>
#   <id> x y z
#   ...
#   $EndNodes
#   $Elements
#   <n_elem>
#   <id> 5 <n_tags> <tag1> <tag2> <n1> ... <n8>
#   ...
#   $EndElements
# ---------------------------------------------------------------------------
function write_msh22(path::String, nodes::Matrix{Float64}, elements::Matrix{Int},
    physical::Vector{Int})
    n_nodes = size(nodes, 1)
    n_elem = size(elements, 1)
    @assert size(nodes, 2) == 3
    @assert size(elements, 2) == 8
    @assert length(physical) == n_elem

    open(path, "w") do io
        println(io, "\$MeshFormat")
        println(io, "2.2 0 8")
        println(io, "\$EndMeshFormat")

        println(io, "\$Nodes")
        println(io, n_nodes)
        for i in 1:n_nodes
            @printf(io, "%d %.17g %.17g %.17g\n", i, nodes[i, 1], nodes[i, 2], nodes[i, 3])
        end
        println(io, "\$EndNodes")

        println(io, "\$Elements")
        println(io, n_elem)
        for e in 1:n_elem
            # 2 etiquetas: física (material) y geométrica (id del elemento).
            @printf(io, "%d 5 2 %d %d", e, physical[e], e)
            for k in 1:8
                @printf(io, " %d", elements[e, k])
            end
            println(io)
        end
        println(io, "\$EndElements")
    end
    println("Escrito: ", path)
end

# ---------------------------------------------------------------------------
# Escritura de un .vtk legacy (unstructured grid) a partir de nodos, elementos
# y etiqueta de material. ParaView lee este formato de forma nativa y fiable
# (a diferencia del lector de .msh, que a veces no muestra nada).
#
# Celdas hexaedro = tipo 12 en VTK.
# ---------------------------------------------------------------------------
function write_vtk(path::String, nodes::Matrix{Float64}, elements::Matrix{Int},
    physical::Vector{Int})
    n_nodes = size(nodes, 1)
    n_elem = size(elements, 1)
    @assert size(nodes, 2) == 3
    @assert size(elements, 2) == 8
    @assert length(physical) == n_elem

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "NedelecHex mesh")
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
        for e in 1:n_elem
            println(io, "12")   # hexaedro
        end
        # Datos escalares por celda: el material (para colorear por εᵣ).
        println(io, "CELL_DATA ", n_elem)
        println(io, "SCALARS material int 1")
        println(io, "LOOKUP_TABLE default")
        for e in 1:n_elem
            println(io, physical[e])
        end
    end
    println("Escrito: ", path)
end

# ---------------------------------------------------------------------------
# 1) Mallas de MATLAB: matlab_nodes_$i.csv + matlab_elems_$i.csv + matlab_eps_$i.csv
# ---------------------------------------------------------------------------
function convert_matlab_meshes()
    for i in 1:5
        nodes_csv = joinpath(OUT_DIR, "matlab_nodes_$i.csv")
        elems_csv = joinpath(OUT_DIR, "matlab_elems_$i.csv")
        eps_csv = joinpath(OUT_DIR, "matlab_eps_$i.csv")
        isfile(nodes_csv) || (println("No existe $nodes_csv, se omite la malla $i"); continue)
        isfile(elems_csv) || (println("No existe $elems_csv, se omite la malla $i"); continue)

        nodes = readdlm(nodes_csv, ',', Float64)
        elems = readdlm(elems_csv, ',', Int)

        # Etiqueta física por elemento: el sólido (material). Si hay eps lo
        # usamos como tag físico para distinguir materiales en Gmsh.
        physical = ones(Int, size(elems, 1))
        if isfile(eps_csv) && filesize(eps_csv) > 0
            eps = readdlm(eps_csv, ',', Float64)
            physical = vec(Int.(round.(eps)))
        end

        write_msh22(joinpath(OUT_DIR, "matlab_mesh_$i.msh"), nodes, elems, physical)
        write_vtk(joinpath(OUT_DIR, "matlab_mesh_$i.vtk"), nodes, elems, physical)
    end
end

# ---------------------------------------------------------------------------
# 2) Mallas generadas con box_mesh (estructuradas, hexaedros alineados).
#    Se exportan varias resoluciones para poder visualizarlas en Gmsh.
# ---------------------------------------------------------------------------
function export_box_meshes()
    # (Lx, Ly, Lz, nx, ny, nz) — cavidad semillena típica del TFG.
    configs = [
        (0.01, 0.001, 0.01, 2, 1, 2),   # 8  elementos
        (0.01, 0.001, 0.01, 4, 1, 4),   # 16 elementos
        (0.01, 0.001, 0.01, 8, 1, 8),   # 64 elementos
    ]
    for (Lx, Ly, Lz, nx, ny, nz) in configs
        mesh = box_mesh(Lx=Lx, Ly=Ly, Lz=Lz, nx=nx, ny=ny, nz=nz)
        n_elem = size(mesh.elements, 1)
        # Etiqueta física: material por defecto (1).
        physical = ones(Int, n_elem)
        base = "box_mesh_$(nx)x$(ny)x$(nz)"
        write_msh22(joinpath(OUT_DIR, base * ".msh"), mesh.nodes, mesh.elements, physical)
        write_vtk(joinpath(OUT_DIR, base * ".vtk"), mesh.nodes, mesh.elements, physical)
    end
end

println("=== Convirtiendo mallas de MATLAB a .msh/.vtk ===")
convert_matlab_meshes()

println("\n=== Exportando mallas generadas con box_mesh ===")
export_box_meshes()

println("\nListo. Todos los .msh y .vtk están en: ", OUT_DIR)
