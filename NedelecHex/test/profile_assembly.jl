#!/usr/bin/env julia
# ===========================================================================
# Perfilador del ensamblaje (dirichlet=false).
# Usa el módulo estándar Profile (sin paquetes extra).
#
# Perfila la malla 3 (4062 dof): representativa de la estructura del bucle
# caliente (mass/stiff elementales, ensamblaje disperso) y ligera, evitando el
# core-dump que provoca perfilar la malla 5 (23124 dof) varias veces.
# dirichlet=false aísla el coste de las matrices elementales (idéntico con o
# sin Dirichlet), que es el objetivo real de optimización.
#
# Salida: output/profile_assembly.txt
# ===========================================================================
using NedelecHex
using Profile
using LinearAlgebra, SparseArrays, Printf

using NedelecHex.Mesh: load_matlab_mesh
using NedelecHex.GlobalDOFs: build_dof_map, dirichlet_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials

BLAS.set_num_threads(1)

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output"))
const OUT_FILE = joinpath(OUTPUT_DIR, "profile_assembly.txt")

# Malla 3 (4062 dof, 800 elem? no: m3 es 4062 dof) — representativa y ligera
mesh = load_matlab_mesh(OUTPUT_DIR, 3)
dof_map = build_dof_map(mesh)
W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

println("Malla 3: n_dof = ", dof_map.n_dofs, ", n_elem = ", size(mesh.elements, 1))

# Función objetivo: el ensamblaje SIN Dirichlet (dirichlet=false). El coste de
# las matrices elementales (mass/stiff) es idéntico con o sin Dirichlet, y
# dirichlet=false aísla el cálculo de las matrices elementales, que es el
# objetivo real de optimización.
assemble() = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=false)

# Calentamiento completo (JIT) ANTES de perfilar.
for _ in 1:6
    assemble()
end

# --- Perfil CPU: bucle de llamadas en régimen estacionario ------------------
println("Perfilando CPU (dirichlet=false, 10 llamadas) ...")
Profile.clear()
@profile for _ in 1:10
    assemble()
end
open(OUT_FILE, "w") do io
    println(io, "=== PERFIL CPU: ensamblaje dirichlet=false, malla 3 (10 llamadas) ===")
    println(io, "n_dof = ", dof_map.n_dofs, ", n_elem = ", size(mesh.elements, 1))
    println(io, "Formato: flat, ordenado por nº de muestras (mincount=10)\n")
    Profile.print(io, format=:flat, sortedby=:count, mincount=10)
end
println("Perfil CPU guardado en ", OUT_FILE)

# --- Perfil de asignaciones (alloc) sobre la misma malla 3 ------------------
println("Perfilando asignaciones (alloc) ...")
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=0.01 assemble()
allocs_report = Profile.Allocs.fetch()
open(OUT_FILE, "a") do io
    println(io, "\n\n=== PERFIL DE ASIGNACIONES (alloc, m3, sample_rate=0.01) ===")
    println(io, "Ordenado por nº de asignaciones (mincount=5)\n")
    Profile.Allocs.print(io, allocs_report; sortedby=:count, mincount=5)
end
println("Perfil de asignaciones guardado en ", OUT_FILE)

println("Listo -> ", OUT_FILE)
