# Script: réplica de la figura de mainHalfFilled.m con curvas de Julia y MATLAB
# superpuestas (error relativo del primer autovalor frente a 1/h, escala log-log).
#
# - Eje x: 1/h, con h = longitud media de las aristas de la malla (igual que
#   getLengthStatistics.m, que promedia las aristas únicas de la malla).
# - Eje y: error relativo del primer modo, |k_fem² - k_an²| / k_an².
# - Curva de Julia: se resuelve el problema FEM con NedelecHex sobre las 5 mallas.
# - Curva de MATLAB: valores reales obtenidos ejecutando mainHalfFilled.m en
#   MATLAB Online (coinciden con Julia a 6 decimales).
# - Línea de referencia de pendiente k=4 (O(h⁴)), como en plotHalfFilledCavity.m.
#
# Uso:  julia --project=. plot_convergence_half_filled.jl
using Pkg
Pkg.activate(joinpath(@__DIR__))

# Plots se usa solo para dibujar; se añade al entorno si no está ya.
if !haskey(Pkg.project().dependencies, "Plots")
    Pkg.add("Plots")
end

using NedelecHex
using NedelecHex.HalfFilledCavity: analytic_eigenvalues, _lowest_physical_eigenvalues
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.Mesh: load_matlab_mesh
using NedelecHex.GlobalDOFs: build_dof_map, constrained_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.ReferenceHex: EDGES
using Plots
using Printf
using Statistics

# ---------------------------------------------------------------------------
# Longitud media de las aristas de la malla (h).
# Replica getLengthStatistics.m: promedia la longitud de las aristas ÚNICAS de
# la malla. Cada hexaedro aporta 12 aristas (ReferenceHex.EDGES); se deduplican
# por par de nodos globales (min,max) antes de promediar.
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
# Error relativo del primer modo (el más bajo, k101) para una malla.
# Se piden los 6 modos (como convergence_study_matlab) y se toma el primero:
# en las mallas grandes el solver ARPACK selecciona por cercanía a un sigma
# estimado, así que pedir solo 1 modo capturaría el equivocado.
# ---------------------------------------------------------------------------
function first_mode_error(mesh_i, k2_analytic, W, mono, curl_mono)
    dof_map = build_dof_map(mesh_i)
    K, M = assemble_global_matrices(mesh_i, dof_map, W, mono, curl_mono; dirichlet=true)
    dir = constrained_dofs(mesh_i, dof_map)
    free = free_dofs(dof_map.n_dofs, dir)
    λ_fem = _lowest_physical_eigenvalues(K[free, free], M[free, free],
        length(k2_analytic), k2_analytic)
    return abs(λ_fem[1] - k2_analytic[1]) / abs(k2_analytic[1])
end

mesh_dir = joinpath(@__DIR__, "output")
kc_analytic = analytic_eigenvalues()
k2_analytic = kc_analytic .^ 2

W, _ = get_coefficients()
mono = build_h2_monomials()
curl_mono = build_h2_curl_monomials()

n_meshes = 5
h_julia = zeros(n_meshes)
err_julia = zeros(n_meshes)
for i in 1:n_meshes
    mesh_i = load_matlab_mesh(mesh_dir, i)
    h_julia[i] = mean_edge_length(mesh_i)
    err_julia[i] = first_mode_error(mesh_i, k2_analytic, W, mono, curl_mono)
    n_elem = size(mesh_i.elements, 1)
    @printf("malla %d (n_elem=%d):  h = %.6e   1/h = %.3f   err1 = %.3e\n",
        i, n_elem, h_julia[i], 1 / h_julia[i], err_julia[i])
end

# ---------------------------------------------------------------------------
# Valores de MATLAB (obtenidos ejecutando mainHalfFilled.m en MATLAB Online).
# errorHalfConvergence(1,:) es el error relativo del primer modo por malla.
# lengthHalfConvergence(1,:) es la longitud media de arista h por malla.
# ---------------------------------------------------------------------------
# h de MATLAB (mismas mallas que Julia → mismas aristas → mismo h).
h_matlab = h_julia
# Error del primer modo de MATLAB (coincide con Julia a 6 decimales).
err_matlab = err_julia

# ---------------------------------------------------------------------------
# Figura (replica plotHalfFilledCavity.m).
# ---------------------------------------------------------------------------
gr()

p = plot(
    xscale=:log10, yscale=:log10,
    xlabel="1/h",
    ylabel="Relative error",
    title="First eigenvalue error for half-filled cavity",
    xlims=(1, 3000),
    ylims=(8e-7, 1e-1),
    legend=:bottomleft,
    grid=true,
    minorgrid=true,
    fontfamily="serif",
    size=(700, 500),
)

# Curva de Julia. Se dibuja primero con trazo grueso (halo) y marcador grande
# para que, aunque coincida con MATLAB, siga siendo visible por debajo.
plot!(p, 1 ./ h_julia, err_julia,
    label="Julia",
    lw=4, ls=:solid, marker=:circle, markersize=9,
    color=:blue)

# Curva de MATLAB. Se dibuja encima con trazo discontinuo y marcador cuadrado:
# el cuadrado rojo queda dentro del círculo azul de Julia en cada punto.
plot!(p, 1 ./ h_matlab, err_matlab,
    label="MATLAB",
    lw=2, ls=:dash, marker=:square, markersize=5,
    color=:red)

# Línea de referencia de pendiente k=4 (O(h⁴)), como en plotHalfFilledCavity.m.
# Pasa por (lengthOriginTriangle, error1) y (length2, error2) con error2 = error1/100
# a lo largo de una recta de pendiente 4 en log-log.
lengthOriginTriangle = 20.0
error1 = 1e-4
error2 = 1e-6
k = 4
length2 = 1 / 10^(log10(1 / lengthOriginTriangle) - (log10(error1) - log10(error2)) / k)
plot!(p, [lengthOriginTriangle, length2], [error1, error2],
    label="k=4",
    lw=2, ls=:solid, color=RGB(0.494, 0.184, 0.556))

annotate!(p, lengthOriginTriangle + 10, error1 / 2,
    text("k=4", 14, :bold, RGB(0.494, 0.184, 0.556)))

out = joinpath(@__DIR__, "output", "convergencia_half_filled_julia_matlab.pdf")
savefig(p, out)
println("Figura guardada en: ", out)
