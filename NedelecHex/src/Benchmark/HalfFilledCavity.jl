module HalfFilledCavity

using LinearAlgebra, SparseArrays, Printf

using ..Mesh: HexMesh, load_hex_mesh
using ..GlobalDOFs: DOFMap, build_dof_map, dirichlet_dofs
using ..Assembly: assemble_global_matrices, free_dofs
using ..CoefficientGenerator: get_coefficients
using ..Monomials: build_h2_monomials, build_h2_curl_monomials

export analytic_kc, analytic_eigenvalues, convergence_study_matlab

# ─── Ecuación trascendente ────────────────────────────────────────────────────
#
# Condición de resonancia para una cavidad rectangular con losa dieléctrica
# (z ∈ [0,h] con εᵣ, z ∈ [h,c] con vacío), pared PEC en todo el contorno.
#
# Derivación: las componentes tangenciales del campo se escriben en cada región
# como sin(kz·z) (vacío en [h,c]) y sin(kz_die·(z-c)) (dieléctrico en [0,h]),
# y la continuidad de Hz en z=h da:
#
#   √(k²−km²)·cot(√(k²−km²)·(c−h)) + √(εᵣk²−km²)·cot(√(εᵣk²−km²)·h) = 0
#
# con km² = (mπ/a)² + (nπ/b)².
#
# Se usa aritmética compleja para manejar de forma uniforme el caso evanescente
# (kz imaginario). Para kz = iκ:
#   kz·cot(kz·L) = (iκ)·cot(iκL) = (iκ)·(−i·coth(κL)) = κ·coth(κL)  (real)
# por lo que tomar la parte real del resultado es siempre correcto.
function _resonance_eq(k, m, n, a, b, c, h, epr)
    k2  = k^2
    km2 = (m * π / a)^2 + (n * π / b)^2

    sq_vac = sqrt(Complex(k2 - km2))
    sq_die = sqrt(Complex(epr * k2 - km2))

    s_vac = sin(sq_vac * (c - h))
    s_die = sin(sq_die * h)

    # Evitar singularidades (cot → ∞) devolviendo Inf
    (abs(s_vac) < 1e-10 || abs(s_die) < 1e-10) && return Inf

    real(sq_vac * cos(sq_vac * (c - h)) / s_vac +
         sq_die * cos(sq_die * h)       / s_die)
end

# ─── Cálculo del valor propio analítico ──────────────────────────────────────
#
# Encuentra la raíz de la ecuación trascendente más cercana a kc_init mediante
# un barrido en malla fina + bisección. La ventana de búsqueda es
# [max(km+ε, kc_init−scan_hw), kc_init+scan_hw].
function analytic_kc(m, n, a, b, c, h, epr, kc_init;
                     scan_hw=300.0, n_scan=30_000, tol=1e-12)
    f   = k -> _resonance_eq(k, m, n, a, b, c, h, epr)
    klo = max(1.0, kc_init - scan_hw)
    khi = kc_init + scan_hw
    ks  = range(klo, khi, n_scan)
    fs  = f.(ks)

    # Recopilar todos los cambios de signo (sin Inf en ninguno de los extremos)
    brackets = Tuple{Float64,Float64}[]
    for i in 1:length(ks)-1
        fi, fj = fs[i], fs[i+1]
        (isfinite(fi) && isfinite(fj) && fi * fj < 0) || continue
        push!(brackets, (ks[i], ks[i+1]))
    end

    isempty(brackets) &&
        error("No se encontró raíz cerca de kc_init=$kc_init para modo ($m,$n)")

    # Elegir la raíz más cercana a kc_init (equivale a vpasolve con punto inicial)
    dists = [abs((lo + hi) / 2 - kc_init) for (lo, hi) in brackets]
    lo, hi = brackets[argmin(dists)]

    for _ in 1:200
        hi - lo < tol && break
        mid = (lo + hi) / 2
        f(lo) * f(mid) <= 0 ? hi = mid : lo = mid
    end
    return (lo + hi) / 2
end

# Calcula los 6 valores propios analíticos del paper para los parámetros dados.
# Los valores iniciales replican los usados en mainHalfFilled.m de MATLAB.
function analytic_eigenvalues(; a=0.01, b=0.001, c=0.01, h=0.005, epr=2.0)
    kw = (m, n, kc0) -> analytic_kc(m, n, a, b, c, h, epr, kc0)
    kc = [kw(1,0,300.0), kw(2,0,540.0), kw(1,0,600.0),
          kw(3,0,750.0), kw(2,0,760.0), kw(1,0,800.0)]
    return kc   # [k101, k201, k102, k301, k202, k103]
end

# ─── εᵣ por elemento ─────────────────────────────────────────────────────────
#
# La losa dieléctrica ocupa z ∈ [0, h_iface]. Un elemento se asigna al material
# si su centroide en z está por debajo de la interfase.
function _element_epsilon(mesh::HexMesh, h_iface::Float64, epr::Float64)
    n_elem  = size(mesh.elements, 1)
    eps_vec = ones(n_elem)
    for e in 1:n_elem
        z_c = sum(mesh.nodes[mesh.elements[e, k], 3] for k in 1:8) / 8
        eps_vec[e] = z_c < h_iface ? epr : 1.0
    end
    return eps_vec
end

# ─── Estudio de convergencia usando mallas exportadas de MATLAB ──────────────
#
# Carga los CSV generados por exportMesh.m y corre el benchmark FEM sobre cada
# malla. Permite comparar directamente con los resultados de mainHalfFilled.m
# usando exactamente la misma discretización del dominio.
function convergence_study_matlab(kc_analytic::Vector{Float64};
                                   mesh_dir::String,
                                   n_meshes::Int=5,
                                   h=0.005, epr=2.0)
    n_modes     = length(kc_analytic)
    k2_analytic = kc_analytic .^ 2
    errors      = zeros(n_modes, n_meshes)

    W, _      = get_coefficients()
    mono      = build_h2_monomials()
    curl_mono = build_h2_curl_monomials()

    sigma = k2_analytic[1] + 1.5 * (kc_analytic[end] - kc_analytic[1])^2

    @printf("\n%-8s  %-8s  %-10s  %s\n", "malla", "n_free", "err_medio", "tiempo")
    @printf("%s\n", "-"^42)

    for i in 1:n_meshes
        nodes_csv = joinpath(mesh_dir, "matlab_nodes_$i.csv")
        elems_csv = joinpath(mesh_dir, "matlab_elems_$i.csv")
        mesh_i  = load_hex_mesh(nodes_csv, elems_csv)
        dof_map = build_dof_map(mesh_i)
        eps_vec = _element_epsilon(mesh_i, h, epr)

        t0 = time()
        K, M = assemble_global_matrices(mesh_i, dof_map, W, mono, curl_mono;
                                        epsilon_r=eps_vec)
        dir  = dirichlet_dofs(mesh_i, dof_map)
        free = free_dofs(dof_map.n_dofs, dir)
        λ_all = eigvals(Symmetric(Matrix(K[free,free])), Symmetric(Matrix(M[free,free])))
        idx   = sortperm(abs.(λ_all .- sigma))[1:n_modes]
        λ_fem = sort(λ_all[idx])
        dt    = time() - t0

        k2_fem      = λ_fem[1:n_modes]
        errors[:,i] = abs.(k2_fem .- k2_analytic) ./ abs.(k2_analytic)
        err_mean    = sum(errors[:,i]) / n_modes
        n_free      = length(free)

        @printf("%-8d  %-8d  %-10.2e  %.1fs\n", i, n_free, err_mean, dt)
    end

    return errors
end

end  # module HalfFilledCavity
