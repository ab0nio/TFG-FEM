using NedelecHex
using Test
using LinearAlgebra, SparseArrays, Printf
# Para ejecutar:  julia --project=. test/runtests.jl

include("test_reference_hex.jl")
include("test_quadrature.jl")
include("test_monomials.jl")
include("test_dof_edge.jl")
include("test_dof_face.jl")
include("test_dof_volume.jl")
include("test_G_matrix.jl")
include("test_jacobian.jl")
include("test_coefficients.jl")
include("test_basis_dual.jl")
include("test_local_matrices.jl")
include("test_single_element.jl")
include("test_two_elements.jl")
include("test_matrices.jl")
include("test_mms.jl")
include("test_analytic_eigenvalues.jl")
include("test_coarse_eigenvalues.jl")
include("test_convergence.jl")
