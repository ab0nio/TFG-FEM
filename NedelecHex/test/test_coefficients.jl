@testset "CoefficientGenerator.jl — Coeficientes W vs Tabla II del paper" begin
    W, G = NedelecHex.CoefficientGenerator.get_coefficients()

    # ── G·W = I₅₄ ─────────────────────────────────────────────────────────────
    @test norm(G * W - I(54)) < 1e-10

    # ── Todas las entradas de W son enteros ──────────────────────────────────
    W_int = round.(Int, W)
    @test maximum(abs.(W .- W_int)) < 1e-8

    # ── Ninguna entrada supera 216 en valor absoluto (Tabla II del paper) ────
    @test maximum(abs.(W_int)) <= 216

    # La comparación exhaustiva columna a columna contra la Tabla II del
    # paper se hace por separado con demo.jl contra la salida de MATLAB
    # (mainGenerateCoefficients.m → tabla_matlab.txt), ya que esos datos no
    # se versionan en el repo. Aquí se comprueban unos valores de referencia
    # puntuales de la Tabla II (filas a1–a2, columnas N₁ y N₂).
    @test W_int[1:2, 1] == [4, -6]
    @test W_int[1:2, 2] == [-2, 6]
end
