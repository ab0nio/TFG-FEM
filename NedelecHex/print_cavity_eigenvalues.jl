# Script temporal: imprime los 6 autovalores FEM (k = sqrt(λ)) por malla de MATLAB
# y los analíticos, para construir la tabla de la memoria.
using Pkg
Pkg.activate(joinpath(@__DIR__))
using NedelecHex
using NedelecHex.HalfFilledCavity: analytic_eigenvalues, _lowest_physical_eigenvalues
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.Mesh: load_matlab_mesh
using NedelecHex.GlobalDOFs: build_dof_map, constrained_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using Printf
using Statistics

mesh_dir = joinpath(@__DIR__, "output")

kc_analytic = analytic_eigenvalues()
k2_analytic = kc_analytic .^ 2
n_modes = length(kc_analytic)

W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

println("=== ANALITICOS ===")
println(join([@sprintf("%.6f", k) for k in kc_analytic], "  "))

for i in 1:5
    mesh_i = load_matlab_mesh(mesh_dir, i)
    dof_map = build_dof_map(mesh_i)
    K, M = assemble_global_matrices(mesh_i, dof_map, W, mono, curl_mono; dirichlet=true)
    dir = constrained_dofs(mesh_i, dof_map)
    free = free_dofs(dof_map.n_dofs, dir)
    λ_fem = _lowest_physical_eigenvalues(K[free, free], M[free, free], n_modes, k2_analytic)
    k_fem = sqrt.(λ_fem)
    n_elem = size(mesh_i.elements, 1)
    n_free = length(free)
    println("=== MALLA $i (n_elem=$n_elem, n_free=$n_free) ===")
    println(join([@sprintf("%.6f", k) for k in k_fem], "  "))
    err = mean(abs.(k_fem .^ 2 .- k2_analytic) ./ k2_analytic)
    println("err_medio = $(@sprintf("%.6e", err))")
end
