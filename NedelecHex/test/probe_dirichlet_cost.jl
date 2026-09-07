#!/usr/bin/env julia
# ===========================================================================
# Sonda: coste real de imponer Dirichlet DENTRO del ensamblaje (dirichlet=true)
# frente a ensamblar la matriz completa y filtrar después (dirichlet=false).
# Aísla el efecto sobre la MISMA malla y en la MISMA ejecución.
# ===========================================================================
using NedelecHex
using BenchmarkTools, LinearAlgebra, SparseArrays, Printf

using NedelecHex.Mesh: load_matlab_mesh
using NedelecHex.GlobalDOFs: build_dof_map, dirichlet_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials

BLAS.set_num_threads(1)

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output"))

# Malla 5 (la más fina, 23124 dof)
mesh = load_matlab_mesh(OUTPUT_DIR, 5)
dof_map = build_dof_map(mesh)
W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

println("Malla 5: n_dof = ", dof_map.n_dofs, ", n_elem = ", size(mesh.elements, 1))

# --- 1. Ensamblaje SIN Dirichlet (matriz completa) ------------------------
assemble_no() = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=false)
assemble_no()  # compila (JIT)
b_no = @benchmark assemble_no() samples=30
t_no = median(b_no).time / 1e6
println("dirichlet=false : ", round(t_no, digits=1), " ms")

# --- 2. Ensamblaje CON Dirichlet dentro (dirichlet=true) -------------------
assemble_yes() = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=true)
assemble_yes()  # compila (JIT)
b_yes = @benchmark assemble_yes() samples=30
t_yes = median(b_yes).time / 1e6
println("dirichlet=true  : ", round(t_yes, digits=1), " ms")

# --- 3. Ensamblaje SIN Dirichlet + filtrado posterior (antigua metodología) -
# El filtrado posterior: extraer el bloque libre-libre (lo que hacía el
# benchmark antiguo para comparar con MATLAB).
K, M = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=false)
free = free_dofs(dof_map.n_dofs, dirichlet_dofs(mesh, dof_map))
filter_post() = (K[free, free], M[free, free])
filter_post()
b_f = @benchmark filter_post() samples=30
t_f = median(b_f).time / 1e6
println("filtrado posterior (K[free,free]) : ", round(t_f, digits=3), " ms")

println()
println("Sobrecoste de dirichlet=true frente a dirichlet=false: ",
    round(100 * (t_yes - t_no) / t_no, digits=1), " %")
println("Sobrecoste del filtrado posterior (antiguo): ",
    round(100 * t_f / t_no, digits=2), " % del ensamblaje sin Dirichlet")
