@testset "GaussQuadrature.jl — Cuadratura de Gauss" begin
    # ── Reglas 1D subyacentes ────────────────────────────────────────────────
    for order in (1, 3, 5)
        w, p = NedelecHex.GaussQuadrature.gauss_segment(order)
        explica("Regla 1D de orden $order: $(length(w)) puntos dentro de [0,1]; Σw = $(sum(w)) — debe ser 1 porque los pesos integran exactamente la función constante 1 (la longitud del intervalo).")
        @test length(w) == length(p)
        @test all(0.0 .<= p .<= 1.0)
        @test isapprox(sum(w), 1.0, atol=1e-12)
    end

    # ── Regla 3D: gauss_hex(4) usa la regla 1D de 3 puntos → 27 puntos ─────────
    weights, points = NedelecHex.GaussQuadrature.gauss_hex(4)
    explica("Regla 3D de orden 4: producto tensorial de la regla 1D de 3 puntos → 3³ = $(length(weights)) puntos; Σw = $(sum(weights)) (volumen del cubo unidad).")
    @test length(weights) == 27
    @test size(points) == (27, 3)

    # Todos los puntos están en [0,1]³
    @test all(0.0 .<= points .<= 1.0)

    # La suma de los pesos es exactamente 1 (partición del volumen unidad)
    @test isapprox(sum(weights), 1.0, atol=1e-12)

    # ── Exactitud: ∫₀¹∫₀¹∫₀¹ ξ²η²ζ² dξdηdζ = (1/3)³ = 1/27 ──────────────────────
    # La regla de 3 puntos por dirección es exacta hasta grado 5, así que
    # integra ξ²η²ζ² (grado 2 en cada variable) sin error.
    f_exact = 1.0 / 27.0
    f_num = sum(weights[k] * points[k,1]^2 * points[k,2]^2 * points[k,3]^2
                for k in eachindex(weights))
    explica("Exactitud polinómica: ∫ξ²η²ζ² dV = (1/3)³ = 1/27 ≈ $(f_exact); la cuadratura da $(f_num). La regla de 3 puntos es exacta hasta grado 5, así que la diferencia es sólo redondeo.")
    @test isapprox(f_num, f_exact, atol=1e-12)

    # ── Regla 2D: gauss_quadrilateral es el producto tensorial de la 1D ────────
    w2, p2 = NedelecHex.GaussQuadrature.gauss_quadrilateral(4)
    explica("Regla 2D (la que usan los DOFs de cara): 3² = $(length(w2)) puntos; Σw = $(sum(w2)) (área del cuadrado unidad).")
    @test length(w2) == 9
    @test isapprox(sum(w2), 1.0, atol=1e-12)
    @test all(0.0 .<= p2 .<= 1.0)
end
