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
