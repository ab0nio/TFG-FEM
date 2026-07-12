@testset "CoefficientGenerator.jl — Matriz G completa" begin
    _, G = NedelecHex.CoefficientGenerator.get_coefficients()

    @test size(G) == (54, 54)

    # ── Unisolvencia ─────────────────────────────────────────────────────────
    # det(G) es formalmente no nulo pero numéricamente diminuto porque las
    # filas mezclan integrales de monomios de grados muy distintos (columnas
    # con escalas muy dispares). El indicador fiable de unisolvencia es que G
    # tiene rango completo y un número de condición moderado en precisión
    # doble (bien por debajo de 1/eps ≈ 10¹⁶).
    @test det(G) != 0.0
    @test rank(G) == 54
    @test cond(G) < 1e5
    explica("Unisolvencia: rank(G) = $(rank(G)) — los 54 DOFs determinan unívocamente cada campo de H₂ — y cond(G) = $(round(cond(G), digits=1)) ≪ 1/eps ≈ 1e16 ⇒ W = G⁻¹ es numéricamente fiable. det(G) solo no es buen indicador: las columnas tienen escalas muy dispares.")

    # ── Estructura de bloques ────────────────────────────────────────────────
    # Cada monomio tiene una única componente cartesiana no nula (columnas
    # 1–18 → ξ̂, 19–36 → η̂, 37–54 → ζ̂). Un DOF sólo puede dar un valor no nulo
    # si su vector de prueba tiene una componente en esa misma dirección.
    col_blocks = (1:18, 19:36, 37:54)

    # Aristas (filas 1–24): cada arista es tangente a un único eje.
    #   ξ̂: aristas 1,3,5,7   η̂: aristas 2,4,6,8   ζ̂: aristas 9,10,11,12
    explica("Ceros estructurales: cada DOF sólo 've' los monomios de su(s) eje(s) compatibles — producto escalar de direcciones ortogonales = 0 exacto (iszero, sin tolerancia) en las 24 filas de arista, 24 de cara y 6 de volumen.")
    edge_axis = Dict(1=>1, 3=>1, 5=>1, 7=>1, 2=>2, 4=>2, 6=>2, 8=>2,
                      9=>3, 10=>3, 11=>3, 12=>3)
    for e in 1:12
        rows = 2*e-1:2*e
        compat = col_blocks[edge_axis[e]]
        incompat = setdiff(1:54, compat)
        @test all(iszero, G[rows, incompat])
    end

    # Caras (filas 25–48): el plano de cada cara determina qué dos de los tres
    # ejes son compatibles (ver Q_FACE en DOFs.jl).
    #   cara 1 (z=0): x,y   cara 2 (z=1): x,y   cara 3 (y=0): x,z
    #   cara 4 (x=1): y,z   cara 5 (y=1): x,z   cara 6 (x=0): y,z
    face_axes = Dict(1=>(1,2), 2=>(1,2), 3=>(1,3), 4=>(2,3), 5=>(1,3), 6=>(2,3))
    for f in 1:6
        base = 24 + 4 * (f - 1)
        rows = base+1:base+4
        compat = union(col_blocks[face_axes[f][1]], col_blocks[face_axes[f][2]])
        incompat = setdiff(1:54, compat)
        @test all(iszero, G[rows, incompat])
    end

    # Volumen (filas 49–54): dos DOFs por eje, en el mismo orden que Q_VOL.
    vol_axis = Dict(1=>1, 2=>1, 3=>2, 4=>2, 5=>3, 6=>3)
    for i in 1:6
        row = 48 + i
        compat = col_blocks[vol_axis[i]]
        incompat = setdiff(1:54, compat)
        @test all(iszero, G[row, incompat])
    end
end
@testset "CoefficientGenerator.jl — Coeficientes W vs Tabla II del paper" begin
    W, G = NedelecHex.CoefficientGenerator.get_coefficients()

    # ── G·W = I₅₄ ─────────────────────────────────────────────────────────────
    res_dual = norm(G * W - I(54))
    explica("Dualidad G·W = I₅₄: residuo de Frobenius = $(res_dual) — cada base Nᵢ cumple gⱼ(Nᵢ) = δᵢⱼ.")
    @test res_dual < 1e-10

    # ── Todas las entradas de W son enteros ──────────────────────────────────
    W_int = round.(Int, W)
    explica("Enteritud: desviación máxima de W respecto a enteros = $(maximum(abs.(W .- W_int))); |coeficiente| máximo = $(maximum(abs.(W_int))) ≤ 216. Que 2916 entradas salgan enteras tras invertir una matriz de integrales en coma flotante es una señal global de corrección: un error en cualquier integral de G rompería la enteritud.")
    @test maximum(abs.(W .- W_int)) < 1e-8

    # ── Ninguna entrada supera 216 en valor absoluto (Tabla II del paper) ────
    @test maximum(abs.(W_int)) <= 216

    # ── Comparación directa contra la Tabla II del paper ─────────────────────
    # test/data/tabla2_paper.txt es la tabla de referencia de coeficientes del
    # paper (exportada de matlab/tableReferenceCoefficients.mat, la que usa
    # mainEvaluateCoefficients.m como "coefficients given in the paper").
    # Columna i = coeficientes de la base Nᵢ sobre los 54 monomios de H₂.
    tabla2 = permutedims(hcat(
        (parse.(Int, split(line)) for line in eachline(joinpath(@__DIR__, "data", "tabla2_paper.txt")))...))
    @test size(tabla2) == (54, 54)
    @test W_int == tabla2

    if TEST_VERBOSE
        println("\n── W calculada (redondeada a entero) ──")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", W_int)
        println("\n\n── Tabla II del paper (referencia) ──")
        show(IOContext(stdout, :limit => true, :displaysize => (30, 140)), "text/plain", tabla2)
        println()
    end
    explica("Comparación directa contra la Tabla II del paper (matlab/tableReferenceCoefficients.mat): entradas distintas = $(count(W_int .!= tabla2)) de $(length(tabla2)).")
end
