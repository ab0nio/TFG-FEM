@testset "BasisEval.jl: Evaluación de las 54 bases vs MATLAB (equivalente a mainEvaluateCoefficients.m)" begin
    # mainEvaluateCoefficients.m evalúa las 54 bases reales en el punto central
    # del hexaedro de referencia (0.5, 0.5, 0.5) de tres formas: (1) simbólico
    # a partir de la Tabla II del paper, (2) el código MATLAB (SystematicBases),
    # (3) Julia: y compara las tres. Aquí reproducimos ese mismo cruce.
    W, mono = W_REF, MONO
    nodes = [NedelecHex.ReferenceHex.NODES[k] for k in 1:8]
    pt = [0.5, 0.5, 0.5]

    N_ref = NedelecHex.BasisEval.eval_bases_ref(W, mono, pt)
    J = NedelecHex.BasisEval.compute_jacobian(nodes, pt)
    N_phys = NedelecHex.BasisEval.piola_transform(N_ref, J)

    explica("Lo que hemos medido antes , el jacobiano en el hexahedro de referencia es J = $(J).
    mainEvaluateCoefficients.m evalúa las 54 bases en el punto central del hexahedro de referencia
    asi que hacemos lo mismo y comparamos")
    @test isapprox(J, I(3), atol=1e-12)

    # Contra la Tabla II simbólica del paper
    # Mismos coeficientes que test_04_coefficients.jl (test/data/tabla2_paper.txt),
    # pero aquí se evalúan las bases en un punto en vez de comparar coeficientes.
    tabla2 = read_matrix(joinpath(@__DIR__, "data", "tabla2_paper.txt"); T=Int)
    N_paper = NedelecHex.BasisEval.eval_bases_ref(tabla2, mono, pt)

    # Contra el código MATLAB (SystematicBases.getBasesEvaluated)
    # test/data/bases_eval_matlab.txt es la salida de mainEvaluateCoefficients.m
    # en el mismo punto (0.5, 0.5, 0.5), copiada de NedelecHex/output/ para que
    # la comparación sea reproducible sin depender de tener MATLAB instalado.
    matlab_eval = read_matrix(joinpath(@__DIR__, "data", "bases_eval_matlab.txt"))
    @test size(matlab_eval) == (54, 3)

    if TEST_VERBOSE
        println("\nBases evaluadas en (0.5, 0.5, 0.5) ")
        println("Tabla II paper (La real):")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", N_paper)
        println("\n\nCódigo MATLAB (SystematicBases):")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", matlab_eval)
        println("\n\nJulia (eval_bases_ref + piola_transform):")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", N_phys)
        println()
    end

    # Comparación triple, igual que mainEvaluateCoefficients.m
    err_paper_matlab = sum(abs.(N_paper .- matlab_eval))
    err_paper_julia = sum(abs.(N_paper .- N_phys))
    err_matlab_julia = sum(abs.(matlab_eval .- N_phys))
    explica("Comparación triple (punto [0.5, 0.5, 0.5]):",
        "\n Error simbólico (Tabla II) vs código MATLAB: $(err_paper_matlab)",
        "\n Error simbólico (Tabla II) vs Julia:         $(err_paper_julia)",
        "\n Error código MATLAB        vs Julia:         $(err_matlab_julia)")

    @test isapprox(N_phys, N_paper, atol=1e-8)
    @test isapprox(N_phys, matlab_eval, atol=1e-8)
    @test isapprox(N_paper, matlab_eval, atol=1e-8)
end
