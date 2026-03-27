using NedelecHex
using Test

@testset "NedelecHex — Tabla II del paper" begin
    passed, W, G = NedelecHex.CoefficientGenerator.run_all_tests(verbose=true)
    @test passed
end