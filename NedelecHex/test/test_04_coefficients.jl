# Prueba CoefficientGenerator.jl. MATLAB: +bases/@SystematicBases/getCoefficients.m;
# la Tabla II de referencia sale de tableReferenceCoefficients.mat (mainEvaluateCoefficients.m:9).
@testset "CoefficientGenerator.jl: Matriz G completa" begin
    G = G_REF

    @test size(G) == (54, 54)

    @test rank(G) == 54
    @test cond(G) < 1e5
    explica("Unisolvencia: rank(G)(número de filas linealmente independientes) = $(rank(G))
    si es 54, los  DOFs determinan unívocamente cada campo de H₂
    cond(G)(número de condicionamiento(mide cuánto la matriz amplifica el error)) = $(round(cond(G), digits=1)) ≪ 1/eps ≈ 1e16 ⇒ W = G⁻¹ es numéricamente fiable.")
end
@testset "CoefficientGenerator.jl: Coeficientes W vs Tabla II del paper" begin
    W, G = W_REF, G_REF

    # G·W = I₅₄
    res_dual = norm(G * W - I(54))
    explica("Dualidad G·W = I₅₄: residuo de Frobenius = $(res_dual) si es pequeño, cada base Nᵢ cumple gⱼ(Nᵢ) = δᵢⱼ
    ya que entonces W es verdaderamente la inversa de G")
    @test res_dual < 1e-10

    # Todas las entradas de W son enteros
    W_int = round.(Int, W)
    explica("Los coeficientes deben ser enteros ,la desviación máxima de W respecto a enteros es $(maximum(abs.(W .- W_int))); el coeficiente máximo en valor absoluto es $(maximum(abs.(W_int))) ≤ 216. Esto no es casualidad: que las 2916 entradas salgan enteras exactas al invertir una matriz de integrales en coma flotante solo pasa si todo el cálculo de G está bien: un error en cualquier integral lo rompería.")
    @test maximum(abs.(W .- W_int)) < 1e-8

    # Ninguna entrada supera 216 en valor absoluto (Tabla II del paper)
    @test maximum(abs.(W_int)) <= 216

    # Comparación directa contra la Tabla II del paper
    # test/data/tabla2_paper.txt es la tabla de referencia de coeficientes del
    # paper (exportada de matlab/tableReferenceCoefficients.mat, la que usa
    # mainEvaluateCoefficients.m como "coefficients given in the paper").
    # Columna i = coeficientes de la base Nᵢ sobre los 54 monomios de H₂.
    tabla2 = read_matrix(joinpath(@__DIR__, "data", "tabla2_paper.txt"); T=Int)
    @test size(tabla2) == (54, 54)
    @test W_int == tabla2

    if TEST_VERBOSE
        println("\nW calculada (redondeada a entero) ")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", W_int)
        println("\n\nTabla II del paper (referencia) ")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", tabla2)
        println()
    end
    explica("Compara contra la Tabla II del paper (matlab/tableReferenceCoefficients.mat): entradas distintas = $(count(W_int .!= tabla2)) de $(length(tabla2)).")
end
