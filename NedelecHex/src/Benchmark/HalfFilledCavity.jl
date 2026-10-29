module HalfFilledCavity
# ymath/plotHalfFilledCavity.m (pero no hace plot solo el error log-log).

using LinearAlgebra, SparseArrays, Printf, DelimitedFiles
using Arpack: eigs

using ..Mesh: HexMesh, load_matlab_mesh
using ..GlobalDOFs: DOFMap, build_dof_map, constrained_dofs
using ..Assembly: assemble_global_matrices, free_dofs
using ..CoefficientGenerator: get_coefficients
using ..Monomials: build_h2_monomials, build_h2_curl_monomials

export analytic_kc, analytic_eigenvalues, convergence_study_matlab

# cavidad rectangular con losa dieléctrica axbxc, altura h, permitividad relativa εᵣ
# (z ∈ [0,h] con εᵣ, z ∈ [h,c] con vacío), pared PEC en todo el contorno.
#
# componentes tangenciales del campo se escriben en cada región
# como sin(kz·z) (vacío en [h,c]) 
# sin(kz_die·(z-c)) (dieléctrico en [0,h]),
# y la continuidad de Hz en z=h da:
#
#   √(k²−km²)·cot(√(k²−km²)·(c−h)) + √(εᵣk²−km²)·cot(√(εᵣk²−km²)·h) = 0
#
# con km² = (mπ/a)² + (nπ/b)².
#
# para manejar de forma uniforme el caso evanescente
# (kz imaginario). Para kz = iκ:
#   kz·cot(kz·L) = (iκ)·cot(iκL) = (iκ)·(−i·coth(κL)) = κ·coth(κL)  (real)
# por lo que tomar la parte real del resultado es válido para ambos casos propagante y evanescente.
function _resonance_eq(k, m, n, a, b, c, h, epr)
    k2 = k^2
    km2 = (m * π / a)^2 + (n * π / b)^2

    sq_vac = sqrt(Complex(k2 - km2))
    sq_die = sqrt(Complex(epr * k2 - km2))

    s_vac = sin(sq_vac * (c - h))
    s_die = sin(sq_die * h)

    # Evitar singularidades (cot → ∞) devolviendo Inf
    (abs(s_vac) < 1e-10 || abs(s_die) < 1e-10) && return Inf

    real(sq_vac * cos(sq_vac * (c - h)) / s_vac +
         sq_die * cos(sq_die * h) / s_die)
end

# Cálculo del valor propio analítico.
# Encuentra la raíz de la ecuación trascendente más cercana a kc_init mediante
# un barrido en malla fina + bisección. La ventana de búsqueda es
# [max(km+ε, kc_init−scan_hw), kc_init+scan_hw].
function analytic_kc(m, n, a, b, c, h, epr, kc_init;
    scan_hw=300.0, n_scan=30_000, tol=1e-12)
    f = k -> _resonance_eq(k, m, n, a, b, c, h, epr)
    klo = max(1.0, kc_init - scan_hw)
    khi = kc_init + scan_hw
    ks = range(klo, khi, n_scan)
    fs = f.(ks)

    # Recopilar todos los cambios de signo (sin Inf en ninguno de los extremos)
    brackets = Tuple{Float64,Float64}[]
    for i in 1:(length(ks)-1)
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
#
# Ojo con el eje: aquí la losa va en z ∈ [0,h], mientras que las mallas de
# meshesHalfFilledTAP2021.mat la tienen cortada perpendicular a x. Comparar
# ambas es válido porque a = c = 0.01 y el corte cae en el plano medio: pasar de
# una a otra es intercambiar x↔z (más una reflexión), una isometría de la
# cavidad, así que el espectro es el mismo. test_10 lo comprueba en las 5 mallas
# en vez de darlo por supuesto.
function analytic_eigenvalues(; a=0.01, b=0.001, c=0.01, h=0.005, epr=2.0)
    kw = (m, n, kc0) -> analytic_kc(m, n, a, b, c, h, epr, kc0)
    kc = [kw(1, 0, 300.0), kw(2, 0, 540.0), kw(1, 0, 600.0),
        kw(3, 0, 750.0), kw(2, 0, 760.0), kw(1, 0, 800.0)]
    return kc   # [k101, k201, k102, k301, k202, k103]
end

# El material ya no se reconstruye aquí: viene etiquetado en la malla
# (mesh.solids / mesh.element_solid) y assemble_global_matrices lo lee solo.
#   - mallas propias: box_mesh + tag_solids_by_plane!(axis=..., cut=...)
# Antes había aquí un `_element_epsilon` que decidía εᵣ por el centroide en z;
# (el azar) y desviaba los autovalores de las mallas gruesas.

# Solver de autovalores.
# Devuelve los n_modes autovalores físicos más bajos de K·x = λ·M·x, dejando
# fuera  λ≈0 (núcleo irrotacional de K).
#
# Mallas pequeñas (n_free ≤ 2000): solver denso + filtro λ > 1.
# Mallas grandes: eigs (ARPACK)  igual que mainHalfFilled.m. Requiere explicittransform=:none
function _lowest_physical_eigenvalues(K, M, n_modes, k2_analytic)
    n = size(K, 1)
    if n <= 2000
        λ_all = eigvals(Symmetric(Matrix(K)), Symmetric(Matrix(M)))
        λ_phys = sort(filter(λ -> λ > 1.0, λ_all))
        return λ_phys[1:n_modes]
    end
    sigma = k2_analytic[1] + 1.5 * (sqrt(k2_analytic[end]) - sqrt(k2_analytic[1]))^2
    λ, _ = eigs(K, M; nev=n_modes, sigma=sigma, which=:LM, explicittransform=:none)
    return sort(real.(λ))
end

# Carga los CSV generados por exportMesh.m y corre el benchmark FEM sobre cada malla
function convergence_study_matlab(kc_analytic::Vector{Float64};
    mesh_dir::String,
    n_meshes::Int=5)
    n_modes = length(kc_analytic)
    k2_analytic = kc_analytic .^ 2
    errors = zeros(n_modes, n_meshes)

    W, _ = get_coefficients()
    mono = build_h2_monomials()
    curl_mono = build_h2_curl_monomials()

    @printf("\n%-8s  %-8s  %-10s  %s\n", "malla", "n_free", "err_medio", "tiempo")
    @printf("%s\n", "-"^42)

    for i in 1:n_meshes
        # decide ni el material ni qué es Dirichlet.
        mesh_i = load_matlab_mesh(mesh_dir, i)
        dof_map = build_dof_map(mesh_i)

        t0 = time()
        # dirichlet=true: impone el PEC por elemento antes de ensamblar, como
        # excluyendo esos DOFs, así que el resultado es idéntico.
        K, M = assemble_global_matrices(mesh_i, dof_map, W, mono, curl_mono;
            dirichlet=true)
        dir = constrained_dofs(mesh_i, dof_map)
        free = free_dofs(dof_map.n_dofs, dir)
        λ_fem = _lowest_physical_eigenvalues(K[free, free], M[free, free],
            n_modes, k2_analytic)
        dt = time() - t0

        k2_fem = λ_fem[1:n_modes]
        errors[:, i] = abs.(k2_fem .- k2_analytic) ./ abs.(k2_analytic)
        err_mean = sum(errors[:, i]) / n_modes
        n_free = length(free)

        @printf("%-8d  %-8d  %-10.2e  %.1fs\n", i, n_free, err_mean, dt)
    end

    return errors
end

end  # module HalfFilledCavity
