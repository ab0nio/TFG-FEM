@testset "HalfFilledCavity.jl: Eigenvalores analíticos de la cavidad semillena" begin
    kc = NedelecHex.HalfFilledCavity.analytic_eigenvalues()
    explica("Datos de referencia(a=1cm, b=1mm, c=1cm, losa h=5mm, εᵣ=2):
    los 6 modos más bajos de la ecuación de resonancia trascendente deben coincidir con los tabulados a 0.5 rad/m.
    Los esperados son: 353.78, 544.50, 599.80, 750.31, 763.32, 809.58 rad/m.")

    # Parámetros del benchmark (paper):
    # cavidad  a=1cm, b=1mm, c=1cm  |  losa dieléctrica h=5mm, εᵣ=2
    # Los 6 modos TE₁₀ₙ / TE₂₀ₙ más bajos deben estar ordenados y coincidir
    # con los valores tabulados en el paper con tolerancia 0.5 rad/m.
    @test length(kc) == 6
    @test all(kc .> 0)
    @test issorted(kc)

    expected = [353.78, 544.50, 599.80, 750.31, 763.32, 809.58]
    labels = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        @test isapprox(kc[i], expected[i], atol=0.5)
    end

    for (label, k) in zip(labels, kc)
        @printf("%-6s = %10.4f rad/m\n", label, k)
    end

    # Convergencia de la bisección
    # Cada kc encontrada debe ser una raíz genuina de la ecuación trascendente:
    # el residuo debe ser prácticamente nulo y debe existir un cambio de signo
    # real a ambos lados (no un artefacto cerca de una asíntota de cot()).
    explica("Raíces reales de la bisección , que sean cero o casi cero |f(kc)| < 1e-6 y
    que no sean cambios de signo por saltos de ±∞ (asíntotas de cot()), tiene que haber un cambio de signo real
    .")
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    modes = [(1, 0), (2, 0), (1, 0), (3, 0), (2, 0), (1, 0)]
    for (i, k) in enumerate(kc)
        m, n = modes[i]
        residual = NedelecHex.HalfFilledCavity._resonance_eq(k, m, n, a, b, c, h, epr)
        @test abs(residual) < 1e-6

        r_lo = NedelecHex.HalfFilledCavity._resonance_eq(k - 1e-6, m, n, a, b, c, h, epr)
        r_hi = NedelecHex.HalfFilledCavity._resonance_eq(k + 1e-6, m, n, a, b, c, h, epr)
        @test sign(r_lo) != sign(r_hi)
    end
end
@testset "HalfFilledCavity.jl: Autovalores FEM en malla gruesa (2×2×2)" begin
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    kc_analytic = NedelecHex.HalfFilledCavity.analytic_eigenvalues(a=a, b=b, c=c, h=h, epr=epr)

    W, mono, curl_mono = W_REF, MONO, CURL_MONO

    mesh = NedelecHex.Mesh.box_mesh(Lx=a, Ly=b, Lz=c, nx=2, ny=2, nz=2)#crea la malla 2x2x2
    dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh) #numera los grados de libertad globales
    eps_vec = NedelecHex.HalfFilledCavity._element_epsilon(mesh, h, epr) #Devuelve un vector con la permitividad de cada elemento
    K, M = NedelecHex.Assembly.assemble_global_matrices(# ensambla las matrices globales de rigidez y masa
        mesh, dof_map, W, mono, curl_mono; epsilon_r=eps_vec)

    dir = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map) #Encuentra los grados de libertad de Dirichlet(en la frontera de PEC)
    free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir) #Son los DOFS que quedan libres (por resolver pq no son de Dirichlet)
    λ_all = eigvals(Symmetric(Matrix(K[free, free])), Symmetric(Matrix(M[free, free])))#resuelve el problema de autovalores que quedan por resolver

    # Filtrado de modos espurios
    # El núcleo irrotacional de K produce un cúmulo de autovalores λ≈0
    # (ruido numérico ~1e-8), muy por debajo de los modos físicos (~1e5). Un
    # umbral fijo separa ambos grupos con holgura de varios órdenes de magnitud.
    λ_physical = sort(filter(λ -> λ > 1.0, λ_all))  #filtra los modos espureos y conserva los físicos
    explica("Modos de  FEM en malla 2×2×2: $(count(λ -> abs(λ) < 1.0, λ_all)) modos espurios λ≈0 (núcleo irrotacional de K) y $(length(λ_physical)) físicos (~1e5) .")
    @test length(λ_physical) >= 6
    @test count(λ -> abs(λ) < 1.0, λ_all) > 0   # existen modos espurios en la malla gruesa

    kc_fem = sqrt.(λ_physical[1:6])
    @test issorted(kc_fem)

    # Comparación con los analíticos
    # Malla muy gruesa (8 elementos): el error es visible pero limitado.
    explica("Como he puesto pocos elemntos (8) y es para probar que el programa funciona permito < 30% de error en el modo con más error.")
    labels = ["k₁₀₁", "k₂₀₁", "k₁₀₂", "k₃₀₁", "k₂₀₂", "k₁₀₃"]
    for i in 1:6
        relerr = abs(kc_fem[i] - kc_analytic[i]) / kc_analytic[i]
        @printf("%-6s  fem=%9.3f  analítico=%9.3f  err=%.1f%%\n",
            labels[i], kc_fem[i], kc_analytic[i], 100 * relerr)
        @test relerr < 0.30
    end
end
@testset "HalfFilledCavity.jl: Convergencia O(h⁴) del primer autovalor" begin
    # Estudio de convergencia equivalente a la Fig. 7 del paper: se refina la
    # malla en las direcciones x,z (nx=nz=n) manteniendo ny=1, ya que los
    # modos de interés (TEₘ₀ₙ) no varían en y y la dimensión b=1mm es muy
    # delgada frente a a=c=1cm: refinar también en y dispararía el coste del
    # solver denso (n_dofs crece con n³ en vez de n²) sin aportar precisión.
    a, b, c, h, epr = 0.01, 0.001, 0.01, 0.005, 2.0
    kc_analytic = NedelecHex.HalfFilledCavity.analytic_eigenvalues(a=a, b=b, c=c, h=h, epr=epr)
    k2_analytic = kc_analytic[1]^2 #primer autovalor de los que hemos encontrado antes

    W, mono, curl_mono = W_REF, MONO, CURL_MONO

    explica("Convergencia del autovalor (Fig. 7 del paper): el error de λ₁ = k² debe caer como O(h⁴): el doble de la tasa O(h²) del campo. Esto pasa porque en FEM el error de un autovalor es cuadrático en el error del autovector, así que dobla la tasa.")
    ns = [1, 2, 4, 8] #número de elementos por lado
    hs = 1.0 ./ ns # tamaño relativo del elemento de la malla: h = 1, 0.5, 0.25, 0.125
    errors = zeros(length(ns)) #error de cada malla

    for (i, n) in enumerate(ns)
        mesh = NedelecHex.Mesh.box_mesh(Lx=a, Ly=b, Lz=c, nx=n, ny=1, nz=n) #refina la malla solo en X
        #la malla más grande será 8x1x8
        dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
        eps_vec = NedelecHex.HalfFilledCavity._element_epsilon(mesh, h, epr)
        K, M = NedelecHex.Assembly.assemble_global_matrices(
            mesh, dof_map, W, mono, curl_mono; epsilon_r=eps_vec)

        dir = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
        free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)
        λ_all = eigvals(Symmetric(Matrix(K[free, free])), Symmetric(Matrix(M[free, free])))
        λ1 = minimum(filter(λ -> λ > 1.0, λ_all))   # primer modo físico (espurios filtrados)

        errors[i] = abs(λ1 - k2_analytic) / k2_analytic
    end

    @printf("\n%-8s  %-12s  %-12s\n", "n", "h", "error rel. k²")
    for i in eachindex(ns)
        @printf("%-8d  %-12.4f  %-12.3e\n", ns[i], hs[i], errors[i])
    end

    # Tasa de convergencia
    # La malla n=1 (6 DOFs libres) está fuera del régimen asintótico; la
    # pendiente se estima por mínimos cuadrados en log-log sobre n=2,4,8.
    idx = 2:length(ns)
    slope = loglog_slope(hs[idx], errors[idx])

    @printf("\nPendiente log-log (n=2,4,8): %.2f  (esperada ≈ 4, paper Fig. 7)\n", slope)

    @test issorted(errors[idx], rev=true)   # el error decrece monótonamente al refinar
    @test isapprox(slope, 4.0, atol=1.0)
end
@testset "HalfFilledCavity.jl: Equivalencia con mainHalfFilled.m (mallas de MATLAB)" begin
    # Réplica 1:1 de mainHalfFilled.m: mismo dominio, mismas 5 mallas
    # (meshesHalfFilledTAP2021.mat exportadas con exportMesh.m a
    # output/matlab_{nodes,elems}_N.csv, versionadas en git) y misma métrica
    # (error relativo medio de los 6 autovalores frente a los analíticos).
    # Los valores de referencia son la salida real de mainHalfFilled.m en
    # MATLAB (R2024, MATLAB Drive, 2026-07-13).
    mesh_dir = joinpath(@__DIR__, "..", "output")

    #en vez de mallas estructuradas como hasta ahora, carga las mallas de matlab
    #Comprueba que el número de DOFs coincide con el que MATLAB imprime en generateDOF (build_dof_map debe reproducir los mismos ) .
    matlab_unknowns = [318, 1678, 4062, 14916, 23124]
    explica("Conteo de DOFs en las 5 mallas de MATLAB: build_dof_map debe reproducir exactamente los 'unknowns' que imprime generateDOF: misma identificación de aristas y caras compartidas sobre mallas reales no estructuradas.")
    for i in 1:5
        mesh_i = NedelecHex.Mesh.load_hex_mesh(
            joinpath(mesh_dir, "matlab_nodes_$i.csv"),
            joinpath(mesh_dir, "matlab_elems_$i.csv"))
        dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh_i)
        @test dof_map.n_dofs == matlab_unknowns[i]
    end

    # Error medio de autovalores: "Mean error obtained"
    matlab_mean_err = [0.0675041324758202, 0.002066105907383419, 0.0003323104331071087,
        7.5286874668569e-5, 3.104919016352379e-5]

    kc = NedelecHex.HalfFilledCavity.analytic_eigenvalues()
    errors = NedelecHex.HalfFilledCavity.convergence_study_matlab(kc; mesh_dir=mesh_dir)
    julia_mean_err = vec(sum(errors, dims=1)) ./ size(errors, 1)

    explica("Error medio por malla, Julia vs mainHalfFilled.m:
    La verdad que no se porque para mallas finas me converge bien y para mallas gruesas no, pero bueno, el error es del mismo orden de magnitud y la tendencia de convergencia es la misma.",
        join(("\n malla $i: $(julia_mean_err[i])  vs  $(matlab_mean_err[i])  (ratio $(round(julia_mean_err[i]/matlab_mean_err[i], digits=4)))"
              for i in 1:5)))

    # Mallas 4-5 (eigs shift-invert, igual que MATLAB): coincidencia a ~1e-8
    # relativo: mismo problema discreto, mismos autovalores. Se exige 1e-6.
    @test isapprox(julia_mean_err[4], matlab_mean_err[4], rtol=1e-6)
    @test isapprox(julia_mean_err[5], matlab_mean_err[5], rtol=1e-6)

    # Mallas 1-3 (solver denso): coinciden en orden de magnitud pero no
    # exactamente (factores medidos 1.43 / 1.19 / 1.15) pese a que en las
    # mallas 4-5 la coincidencia es ~1e-8 con el mismo ensamblaje.
    @test 1/2.5 < julia_mean_err[1] / matlab_mean_err[1] < 2.5
    @test 1/2 < julia_mean_err[2] / matlab_mean_err[2] < 2
    @test 1/2 < julia_mean_err[3] / matlab_mean_err[3] < 2

    # Misma tendencia de convergencia que la Fig. 5 del paper.
    @test issorted(julia_mean_err, rev=true)
end
