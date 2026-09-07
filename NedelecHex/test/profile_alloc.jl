#!/usr/bin/env julia
# ===========================================================================
# Perfilador de ASIGNACIONES (alloc) del ensamblaje (dirichlet=true).
# Versión mínima y rápida: malla 1 (más pequeña), 1 sola llamada, para evitar
# el core-dump / lentitud extrema de perfilar mallas grandes con sample_rate
# alto. Escribe la salida con flush para que sobreviva a interrupciones.
#
# Salida: output/profile_alloc.txt
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
const OUT_FILE = joinpath(OUTPUT_DIR, "profile_alloc.txt")

mesh = load_matlab_mesh(OUTPUT_DIR, 1)
dof_map = build_dof_map(mesh)
W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

println("Malla 1: n_dof = ", dof_map.n_dofs, ", n_elem = ", size(mesh.elements, 1))
flush(stdout)

assemble() = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=true)

# Calentamiento completo (JIT) ANTES de perfilar.
for _ in 1:6
    assemble()
end
println("Calentamiento listo. Perfilando alloc (m1, 1 llamada, sample_rate=0.1) ...")
flush(stdout)

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=0.1 assemble()
allocs_report = Profile.Allocs.fetch()
open(OUT_FILE, "w") do io
    println(io, "=== PERFIL DE ASIGNACIONES (alloc, m1, 1 llamada, sample_rate=0.1) ===")
    println(io, "n_dof = ", dof_map.n_dofs, ", n_elem = ", size(mesh.elements, 1))
    println(io, "Ordenado por nº de asignaciones (mincount=5)\n")
    Profile.Allocs.print(io, allocs_report; sortedby=:count, mincount=5)
end
println("Perfil de asignaciones guardado en ", OUT_FILE)
flush(stdout)
println("Listo -> ", OUT_FILE)
