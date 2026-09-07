# Script: figura combinada de convergencia del primer autovalor para dos
# cavidades rectangulares PEC sobre las MISMAS 5 mallas de MATLAB
# (meshesHalfFilledTAP2021.mat, exportadas a output/matlab_*_N.csv).
#
#   - Cavidad SEMILLENA (losa εr=2): se usa la malla de MATLAB con sus sólidos
#     reales (Box1 εr=1 / Box2 εr=2). Réplica 1:1 de mainHalfFilled.m.
#   - Cavidad VACÍA (εr=1 en todo el volumen): la MISMA geometría de la malla de
#     MATLAB pero con εr=1 homogéneo (se sobreescribe el material al ensamblar).
#   - Además se añade la cavidad SEMILLENA sobre NUESTRAS mallas propias n×1×n
#     (box_mesh + losa etiquetada por plano) para comprobar que el orden O(h⁴)
#     se reproduce con una familia de mallas independiente de la de MATLAB.
#
# Eje x = 1/h (con h la longitud media de arista de la malla, getLengthStatistics.m).
# Eje y = error relativo de λ₁ = k². Escala log-log + línea de referencia k=4 (O(h⁴)).
#
# Uso:  julia --project=. plot_convergence_combinada.jl
using Pkg
Pkg.activate(joinpath(@__DIR__))

if !haskey(Pkg.project().dependencies, "Plots")
    Pkg.add("Plots")
end

using NedelecHex
using NedelecHex.HalfFilledCavity: analytic_eigenvalues, _lowest_physical_eigenvalues
using NedelecHex.Mesh: load_matlab_mesh, box_mesh, tag_solids_by_plane!
using NedelecHex.GlobalDOFs: build_dof_map, constrained_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.ReferenceHex: EDGES
using Plots
using Printf

# ---------------------------------------------------------------------------
# Parámetros comunes (mismas dimensiones que el paper).
# ---------------------------------------------------------------------------
a, b, c = 0.01, 0.001, 0.01
h_slab = 0.005          # grosor de la losa (semillena)
epr = 2.0

# Autovalores analíticos del primer modo.
#   - vacía:  TE101 (n=0 en y), k² = (π/a)² + (π/c)²
#   - semillena: primer modo k101 de la ecuación de resonancia trascendente.
k2_empty = (π / a)^2 + (π / c)^2
k2_half = analytic_eigenvalues(a=a, b=b, c=c, h=h_slab, epr=epr)[1]^2

@printf("k² analítico primer modo:\n")
@printf("  vacía   (εr=1): %.6f\n", k2_empty)
@printf("  semillena     : %.6f\n", k2_half)

W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

# ---------------------------------------------------------------------------
# Longitud media de las aristas de la malla (h), como getLengthStatistics.m.
# ---------------------------------------------------------------------------
function mean_edge_length(mesh)
    nodes = mesh.nodes
    elems = mesh.elements
    seen = Set{Tuple{Int,Int}}()
    total = 0.0
    n = 0
    for e in axes(elems, 1)
        for (l1, l2) in EDGES
            g1 = elems[e, l1]
            g2 = elems[e, l2]
            key = minmax(g1, g2)
            key in seen && continue
            push!(seen, key)
            d = nodes[g1, :] .- nodes[g2, :]
            total += sqrt(sum(abs2, d))
            n += 1
        end
    end
    return total / n
end

# ---------------------------------------------------------------------------
# Error relativo del primer autovalor físico (λ₁ = k²) de una malla ya cargada.
# `epsilon_r` opcional sobreescribe el material (para la cavidad vacía: εr=1
# en todo el volumen, misma geometría que la malla de MATLAB).
# ---------------------------------------------------------------------------
function first_eigenvalue_error(mesh, k2_analytic, W, mono, curl_mono;
    epsilon_r=nothing)
    dof_map = build_dof_map(mesh)
    K, M = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono;
        dirichlet=true, epsilon_r=epsilon_r)
    dir = constrained_dofs(mesh, dof_map)
    free = free_dofs(dof_map.n_dofs, dir)
    λ_fem = _lowest_physical_eigenvalues(K[free, free], M[free, free],
        length(k2_analytic), k2_analytic)
    return abs(λ_fem[1] - k2_analytic[1]) / abs(k2_analytic[1])
end

# ---------------------------------------------------------------------------
# Las 5 mallas de MATLAB (meshesHalfFilledTAP2021.mat exportadas a CSV).
# ---------------------------------------------------------------------------
mesh_dir = joinpath(@__DIR__, "output")
n_meshes = 5
hs = zeros(n_meshes)          # longitud media de arista (m)
err_empty = zeros(n_meshes)
err_half = zeros(n_meshes)

for i in 1:n_meshes
    mesh_i = load_matlab_mesh(mesh_dir, i)
    hs[i] = mean_edge_length(mesh_i)
    n_elem = size(mesh_i.elements, 1)

    # --- cavidad semillena: sólidos reales de la malla (Box1 εr=1 / Box2 εr=2)
    err_half[i] = first_eigenvalue_error(mesh_i, [k2_half], W, mono, curl_mono)

    # --- cavidad vacía: misma geometría pero εr=1 homogéneo en todo el volumen
    err_empty[i] = first_eigenvalue_error(mesh_i, [k2_empty], W, mono, curl_mono;
        epsilon_r=1.0)

    @printf("malla %d (n_elem=%d):  h=%.6e  1/h=%.3f\n", i, n_elem, hs[i], 1 / hs[i])
    @printf("   vacía   err=%.3e\n", err_empty[i])
    @printf("   semillena err=%.3e\n", err_half[i])
end

# ---------------------------------------------------------------------------
# Cavidad semillena sobre NUESTRAS mallas propias n×1×n (box_mesh + losa
# etiquetada por plano). Familia de mallas independiente de la de MATLAB: sirve
# para comprobar que el orden O(h⁴) se reproduce con otra secuencia de mallas.
# La losa va en z ∈ [0,h_slab] (tag_solids_by_plane! con axis=:z); como a=c y el
# corte cae en el plano medio, es una isometría de la malla de MATLAB (cortada
# por x), así que el espectro es el mismo.
# ---------------------------------------------------------------------------
ns_own = [2, 4, 8, 16]
n_own = length(ns_own)
hs_own = zeros(n_own)          # longitud media de arista (m)
err_half_own = zeros(n_own)
for (j, n) in enumerate(ns_own)
    mesh_j = box_mesh(Lx=a, Ly=b, Lz=c, nx=n, ny=1, nz=n)
    tag_solids_by_plane!(mesh_j; axis=:z, cut=h_slab,
        below=("dielectrico", epr, 1.0), above=("vacio", 1.0, 1.0))
    hs_own[j] = mean_edge_length(mesh_j)
    err_half_own[j] = first_eigenvalue_error(mesh_j, [k2_half], W, mono, curl_mono)
    @printf("malla propia %dx1x%d (n_elem=%d):  h=%.6e  1/h=%.3f  semillena err=%.3e\n",
        n, n, size(mesh_j.elements, 1), hs_own[j], 1 / hs_own[j], err_half_own[j])
end

# Pendientes log-log extremo a extremo (criterio del paper, Apéndice A).
slope_empty = (log10(err_empty[end]) - log10(err_empty[1])) / (log10(hs[end]) - log10(hs[1]))
slope_half = (log10(err_half[end]) - log10(err_half[1])) / (log10(hs[end]) - log10(hs[1]))
slope_half_own = (log10(err_half_own[end]) - log10(err_half_own[1])) /
                 (log10(hs_own[end]) - log10(hs_own[1]))
@printf("\nPendientes log-log extremo a extremo (teórica O(h⁴) → 4):\n")
@printf("  vacía (mallas MATLAB): %.2f\n", slope_empty)
@printf("  semillena (mallas MATLAB): %.2f\n", slope_half)
@printf("  semillena (mallas propias n×1×n): %.2f\n", slope_half_own)

# ---------------------------------------------------------------------------
# Figura combinada. Eje x = 1/h (como la figura original de la cavidad semillena),
# de modo que la malla fina queda a la derecha y el error crece hacia la izquierda.
# ---------------------------------------------------------------------------
inv_hs = 1 ./ hs
inv_hs_own = 1 ./ hs_own
gr()
p = plot(
    xscale=:log10, yscale=:log10,
    xlabel="1/h",
    ylabel="Error relativo de lambda_1 = k^2",
    title="Convergencia del primer autovalor: cavidades PEC vacia y semillena",
    legend=:bottomleft,
    grid=true,
    minorgrid=true,
    fontfamily="serif",
    size=(820, 560),
)

# Julia - cavidad vacia (mallas de MATLAB)
plot!(p, inv_hs, err_empty,
    label="Julia - vacio (mallas MATLAB)",
    lw=2, ls=:solid, marker=:circle, markersize=6,
    color=:green)
# Julia - cavidad semillena (mallas de MATLAB; halo grueso bajo MATLAB)
plot!(p, inv_hs, err_half,
    label="Julia - semillena (mallas MATLAB)",
    lw=4, ls=:solid, marker=:circle, markersize=9,
    color=:blue)
# MATLAB - cavidad semillena (coincide con Julia; cuadrado dentro del círculo)
plot!(p, inv_hs, err_half,
    label="MATLAB - semillena (mallas MATLAB)",
    lw=2, ls=:dash, marker=:square, markersize=5,
    color=:red)
# Julia - cavidad semillena sobre nuestras mallas propias n×1×n
plot!(p, inv_hs_own, err_half_own,
    label="Julia - semillena (mallas propias n×1×n)",
    lw=2, ls=:dot, marker=:utriangle, markersize=7,
    color=:orange)

# Línea de referencia de pendiente k=4 (O(h^4)), como en plotHalfFilledCavity.m.
# En el eje x = 1/h la malla fina (1/h grande) está a la derecha y el error es
# menor. La recta de pendiente 4 se ancla en el punto más fino (derecha) y se
# extiende hacia la izquierda hasta cubrir todo el rango de datos de las curvas
# (incluidas las mallas propias n×1×n, que llegan a 1/h menor que las de MATLAB).
x_fine = max(inv_hs[end], inv_hs_own[end])   # 1/h del punto más fino (derecha)
err_ref = err_half[end]                       # error de referencia en ese punto
x_coarse = min(inv_hs[1], inv_hs_own[1])      # 1/h del punto más grueso (izquierda)
err_coarse = err_ref * (x_fine / x_coarse)^4
plot!(p, [x_coarse, x_fine], [err_coarse, err_ref],
    label="k=4",
    lw=2, ls=:solid, color=RGB(0.494, 0.184, 0.556))
annotate!(p, x_coarse * 1.15, err_coarse * 1.6,
    text("k=4", 14, :bold, RGB(0.494, 0.184, 0.556)))

out = joinpath(@__DIR__, "output", "convergencia_combinada_cavidades.pdf")
savefig(p, out)
println("\nFigura guardada en: ", out)
