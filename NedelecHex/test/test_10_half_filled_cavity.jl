@testset "HalfFilledCavity.jl — Eigenvalores analíticos de la cavidad semillena" begin
    kc = NedelecHex.HalfFilledCavity.analytic_eigenvalues()
    explica("Referencia analítica del benchmark (a=1cm, b=1mm, c=1cm, losa h=5mm, εᵣ=2): los 6 modos más bajos de la ecuación de resonancia trascendente deben coincidir con los tabulados a 0.5 rad/m.")

    # Parámetros del benchmark (paper Amor-Martin 2023):
    #   cavidad  a=1cm, b=1mm, c=1cm  |  losa dieléctrica h=5mm, εᵣ=2
    # Los 6 modos TE₁₀ₙ / TE₂₀ₙ más bajos deben estar ordenados y coincidir
    # con los valores tabulados en el paper con tolerancia 0.5 rad/m.
    @test length(kc) == 6
    @test all(kc .> 0)
    @test issorted(kc)

    expected = [353.78, 544.50, 599.80, 750.31, 763.32, 809.58]
    labels   = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        @test isapprox(kc[i], expected[i], atol=0.5)
    end

    for (label, k) in zip(labels, kc)
        @printf("%-6s = %10.4f rad/m\n", label, k)
    end

    # ── Convergencia de la bisección ─────────────────────────────────────────
    # Cada kc encontrada debe ser una raíz genuina de la ecuación trascendente:
    # el residuo debe ser prácticamente nulo y debe existir un cambio de signo
    # real a ambos lados (no un artefacto cerca de una asíntota de cot()).
    explica("Raíces genuinas de la bisección: residuo |f(kc)| < 1e-6 y cambio de signo real a ambos lados — descarta que la 'raíz' sea una asíntota de cot(), donde también hay cruce aparente.")
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    modes = [(1,0), (2,0), (1,0), (3,0), (2,0), (1,0)]
    for (i, k) in enumerate(kc)
        m, n = modes[i]
        residual = NedelecHex.HalfFilledCavity._resonance_eq(k, m, n, a, b, c, h, epr)
        @test abs(residual) < 1e-6

        r_lo = NedelecHex.HalfFilledCavity._resonance_eq(k - 1e-6, m, n, a, b, c, h, epr)
        r_hi = NedelecHex.HalfFilledCavity._resonance_eq(k + 1e-6, m, n, a, b, c, h, epr)
        @test sign(r_lo) != sign(r_hi)
    end
end
@testset "HalfFilledCavity.jl — Autovalores FEM en malla gruesa (2×2×2)" begin
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    kc_analytic = NedelecHex.HalfFilledCavity.analytic_eigenvalues(a=a, b=b, c=c, h=h, epr=epr)

    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    mesh    = NedelecHex.Mesh.box_mesh(Lx=a, Ly=b, Lz=c, nx=2, ny=2, nz=2)
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
    eps_vec = NedelecHex.HalfFilledCavity._element_epsilon(mesh, h, epr)
    K, M    = NedelecHex.Assembly.assemble_global_matrices(
                  mesh, dof_map, W, mono, curl_mono; epsilon_r=eps_vec)

    dir  = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
    free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)
    λ_all = eigvals(Symmetric(Matrix(K[free, free])), Symmetric(Matrix(M[free, free])))

    # ── Filtrado de modos espurios ────────────────────────────────────────────
    # El núcleo irrotacional de K produce un cúmulo de autovalores λ≈0
    # (ruido numérico ~1e-8), muy por debajo de los modos físicos (~1e5). Un
    # umbral fijo separa ambos grupos con holgura de varios órdenes de magnitud.
    λ_physical = sort(filter(λ -> λ > 1.0, λ_all))
    explica("Espectro FEM en malla 2×2×2: $(count(λ -> abs(λ) < 1.0, λ_all)) modos espurios λ≈0 (núcleo irrotacional de K) y $(length(λ_physical)) físicos (~1e5); el umbral λ>1 los separa con varios órdenes de magnitud de holgura.")
    @test length(λ_physical) >= 6
    @test count(λ -> abs(λ) < 1.0, λ_all) > 0   # existen modos espurios en la malla gruesa

    kc_fem = sqrt.(λ_physical[1:6])
    @test issorted(kc_fem)

    # ── Comparación con los analíticos ────────────────────────────────────────
    # Malla muy gruesa (8 elementos): el error es visible pero limitado.
    explica("Con sólo 8 elementos se exige error < 30% por modo: es una verificación de plausibilidad del pipeline completo; la precisión real se mide en el estudio de convergencia siguiente.")
    labels = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        relerr = abs(kc_fem[i] - kc_analytic[i]) / kc_analytic[i]
        @printf("%-6s  fem=%9.3f  analítico=%9.3f  err=%.1f%%\n",
                labels[i], kc_fem[i], kc_analytic[i], 100 * relerr)
        @test relerr < 0.30
    end
end
@testset "HalfFilledCavity.jl — Convergencia O(h⁴) del primer autovalor" begin
    # Estudio de convergencia equivalente a la Fig. 7 del paper: se refina la
    # malla en las direcciones x,z (nx=nz=n) manteniendo ny=1, ya que los
    # modos de interés (TEₘ₀ₙ) no varían en y y la dimensión b=1mm es muy
    # delgada frente a a=c=1cm — refinar también en y dispararía el coste del
    # solver denso (n_dofs crece con n³ en vez de n²) sin aportar precisión.
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    kc_analytic = NedelecHex.HalfFilledCavity.analytic_eigenvalues(a=a, b=b, c=c, h=h, epr=epr)
    k2_analytic = kc_analytic[1]^2

    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    explica("Convergencia del autovalor (Fig. 7 del paper): el error de λ₁ = k² debe caer como O(h⁴) — el doble de la tasa O(h²) del campo, propiedad estándar de los autovalores en FEM (su error es cuadrático en el error del autovector).")
    ns     = [1, 2, 4, 8]
    hs     = 1.0 ./ ns
    errors = zeros(length(ns))

    for (i, n) in enumerate(ns)
        mesh    = NedelecHex.Mesh.box_mesh(Lx=a, Ly=b, Lz=c, nx=n, ny=1, nz=n)
        dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
        eps_vec = NedelecHex.HalfFilledCavity._element_epsilon(mesh, h, epr)
        K, M    = NedelecHex.Assembly.assemble_global_matrices(
                      mesh, dof_map, W, mono, curl_mono; epsilon_r=eps_vec)

        dir   = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
        free  = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)
        λ_all = eigvals(Symmetric(Matrix(K[free, free])), Symmetric(Matrix(M[free, free])))
        λ1    = minimum(filter(λ -> λ > 1.0, λ_all))   # primer modo físico (espurios filtrados)

        errors[i] = abs(λ1 - k2_analytic) / k2_analytic
    end

    @printf("\n%-8s  %-12s  %-12s\n", "n", "h", "error rel. k²")
    for i in eachindex(ns)
        @printf("%-8d  %-12.4f  %-12.3e\n", ns[i], hs[i], errors[i])
    end

    # ── Tasa de convergencia ──────────────────────────────────────────────────
    # La malla n=1 (6 DOFs libres) está fuera del régimen asintótico; la
    # pendiente se estima por mínimos cuadrados en log-log sobre n=2,4,8.
    idx = 2:length(ns)
    x   = log.(hs[idx])
    y   = log.(errors[idx])
    x̄, ȳ = sum(x)/length(x), sum(y)/length(y)
    slope = sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄).^2)

    @printf("\nPendiente log-log (n=2,4,8): %.2f  (esperada ≈ 4, paper Fig. 7)\n", slope)

    @test issorted(errors[idx], rev=true)   # el error decrece monótonamente al refinar
    @test isapprox(slope, 4.0, atol=1.0)
end
