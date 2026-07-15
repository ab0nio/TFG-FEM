using NedelecHex
using Test
using LinearAlgebra, SparseArrays, Printf
# Para ejecutar:  julia --project=. test/runtests.jl
# Modo verboso:   NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl
# Cada bloque imprime QUÉ comprueba, POR QUÉ el valor esperado
# es el correcto y las magnitudes medidas, además del pasa/falla.
const TEST_VERBOSE = lowercase(get(ENV, "NEDELEC_TEST_VERBOSE", "")) in ("1", "true", "yes")
explica(args...) = (TEST_VERBOSE && println("  ▸ ", args..., "\n"); nothing)

# Fixtures compartidas por toda la suite
# W/G y los monomios son inmutables y G son 54×54 integrales relativamente
# caras: se calculan una sola vez aquí en lugar de en cada testset. Los tests
# 02, 04, 05 y 07 validan estos mismos objetos contra el paper y MATLAB.
const W_REF, G_REF = NedelecHex.CoefficientGenerator.get_coefficients()
const MONO         = NedelecHex.Monomials.build_h2_monomials()
const CURL_MONO    = NedelecHex.Monomials.build_h2_curl_monomials()

# Matriz desde texto plano (una fila por línea, columnas separadas por
# espacios), para los ficheros de referencia de test/data/.
read_matrix(path; T=Float64) =
    permutedims(hcat((parse.(T, split(line)) for line in eachline(path))...))

# Pendiente de log(err) vs log(h) por mínimos cuadrados: el orden de
# convergencia observado en un estudio de refinamiento.
function loglog_slope(hs, errs)
    x, y = log.(hs), log.(errs)
    x̄, ȳ = sum(x)/length(x), sum(y)/length(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄).^2)
end
#
# Orden lógico: cada nivel depende solo de los anteriores.
# 01-02  infraestructura matemática (cuadratura, espacio H₂)
# 03-04  funcionales DOF y coeficientes (Tabla II del paper)
# 05     elemento de referencia (Jacobiano, base dual, Kᵉ/Mᵉ)
# 06     malla (permutación de nodos de Gmsh)
# 07     bases evaluadas vs MATLAB (equivalente a mainEvaluateCoefficients.m)
# 08     ensamblaje global en progresión 1 → 2 → 8 elementos
# 10     verificación física (cavidad semillena)

for f in (
    "test_01_quadrature.jl",
    "test_02_monomials.jl",
    "test_03_dofs.jl",
    "test_04_coefficients.jl",
    "test_05_reference_element.jl",
    "test_06_gmsh_node_order.jl",
    "test_07_bases_eval.jl",
    "test_08_assembly.jl",
    "test_10_half_filled_cavity.jl",
)
    TEST_VERBOSE && println("\n", "═"^70, "\n▶ Ejecutando ", f, "\n", "═"^70)
    include(f)
end
