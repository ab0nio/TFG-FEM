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
    @test length(λ_physical) >= 6
    @test count(λ -> abs(λ) < 1.0, λ_all) > 0   # existen modos espurios en la malla gruesa

    kc_fem = sqrt.(λ_physical[1:6])
    @test issorted(kc_fem)

    # ── Comparación con los analíticos ────────────────────────────────────────
    # Malla muy gruesa (8 elementos): el error es visible pero limitado.
    labels = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        relerr = abs(kc_fem[i] - kc_analytic[i]) / kc_analytic[i]
        @printf("%-6s  fem=%9.3f  analítico=%9.3f  err=%.1f%%\n",
                labels[i], kc_fem[i], kc_analytic[i], 100 * relerr)
        @test relerr < 0.30
    end
end
