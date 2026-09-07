# Script: estudio de convergencia de la cavidad rectangular PEC VACÍA (homogénea).
#
# Cavidad a×b×c con paredes PEC y sin dieléctrico. Los autovalores son exactos y
# triviales (separables):
#     k² = (mπ/a)² + (nπ/b)² + (pπ/c)²
# El modo más bajo (b muy fino → n=0) es el TE101: k² = (π/a)² + (π/c)².
#
# Se refina la malla en x,z (nx=nz=n, ny=1) y se mide el error relativo de
# λ₁ = k² frente a h. La pendiente log-log debe ser O(h⁴) (elemento de Nédélec
# de orden 2). La cavidad vacía converge muy rápido, así que se descarta la
# malla 1×1×1 (6 DOFs, fuera del régimen asintótico) y se refina hasta n=16.
#
# Uso:  julia --project=. rect_cavity_convergence.jl
using Pkg
Pkg.activate(joinpath(@__DIR__))

if !haskey(Pkg.project().dependencies, "Plots")
    Pkg.add("Plots")
end

using NedelecHex
using NedelecHex.Mesh: box_mesh
using NedelecHex.GlobalDOFs: build_dof_map, dirichlet_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials
using NedelecHex.CoefficientGenerator: get_coefficients
using Plots
using Printf
using LinearAlgebra
using Statistics

# ---------------------------------------------------------------------------
# Parámetros de la cavidad (mismas dimensiones que la semillena, pero vacía).
# ---------------------------------------------------------------------------
a, b, c = 0.01, 0.001, 0.01

# Autovalor analítico del primer modo TE101 (n=0 en la dirección fina y).
k2_analytic = (π / a)^2 + (π / c)^2
k_analytic = sqrt(k2_analytic)
@printf("Cavidad vacía %.3g × %.3g × %.3g m\n", a, b, c)
@printf("Primer modo TE101: k = %.6f rad/m, k² = %.6f\n", k_analytic, k2_analytic)

W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

# ---------------------------------------------------------------------------
# Estudio de convergencia n×1×n.
# ---------------------------------------------------------------------------
ns = [2, 4, 8, 16]
hs = 1.0 ./ ns          # tamaño relativo del elemento (h = 1/n)
errors = zeros(length(ns))
for (i, n) in enumerate(ns)
    mesh = box_mesh(Lx=a, Ly=b, Lz=c, nx=n, ny=1, nz=n)
    dof_map = build_dof_map(mesh)
    K, M = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono)
    dir = dirichlet_dofs(mesh, dof_map)
    free = free_dofs(dof_map.n_dofs, dir)
    n_free = length(free)
    # Resolución densa en todas las mallas (hasta n=16, 3906 DOFs). La cavidad
    # vacía tiene un espectro denso cerca de k² y un gran núcleo irrotacional en
    # λ≈0, lo que hace poco fiable a eigs/ARPACK con shift; el solver denso con
    # filtro λ>1 es robusto y factible hasta ~5000 DOFs.
    λ_all = eigvals(Symmetric(Matrix(K[free, free])), Symmetric(Matrix(M[free, free])))
    λ1 = minimum(filter(λ -> λ > 1.0, λ_all))   # primer modo físico (espurios filtrados)
    errors[i] = abs(λ1 - k2_analytic) / k2_analytic
    @printf("malla %dx1x%d:  h=%.3f  n_free=%d  λ1=%.6f  err=%.3e\n",
        n, n, hs[i], n_free, λ1, errors[i])
end

# Pendiente log-log extremo a extremo (criterio del paper, Apéndice A).
slope = (log10(errors[end]) - log10(errors[1])) / (log10(hs[end]) - log10(hs[1]))
@printf("\nPendiente log-log extremo a extremo: %.2f  (teórica O(h⁴) → 4)\n", slope)

# ---------------------------------------------------------------------------
# Figura: error relativo de λ₁ frente a h (log-log) + línea de referencia O(h⁴).
# ---------------------------------------------------------------------------
gr()
p = plot(
    xscale=:log10, yscale=:log10,
    xlabel="h",
    ylabel="Relative error of λ₁ = k²",
    title="Empty rectangular PEC cavity — first eigenvalue convergence",
    legend=:bottomright,
    grid=true,
    minorgrid=true,
    fontfamily="serif",
    size=(700, 500),
)
plot!(p, hs, errors,
    label="Julia (empty cavity)",
    lw=2, ls=:solid, marker=:circle, markersize=6,
    color=:blue)

# Línea de referencia de pendiente 4 (O(h⁴)) anclada en el punto más fino.
h_ref = hs[end]
err_ref = errors[end]
# recta de pendiente 4 en log-log: err = err_ref * (h/h_ref)^4
h_lo, h_hi = hs[1], hs[end]
plot!(p, [h_lo, h_hi], [err_ref * (h_lo / h_ref)^4, err_ref],
    label="O(h⁴)",
    lw=2, ls=:dash, color=RGB(0.494, 0.184, 0.556))

out = joinpath(@__DIR__, "output", "convergencia_cavidad_vacia.pdf")
savefig(p, out)
println("Figura guardada en: ", out)
