# Equivalente en MATLAB: no hay suite de tests; los drivers son
# mainGenerateCoefficients.m, mainEvaluateCoefficients.m y mainHalfFilled.m.
using NedelecHex
using Test
using LinearAlgebra, SparseArrays, Printf
# Para ejecutar:  julia --project=. test/runtests.jl
# Modo verboso:   NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl
# Omitir tests:   NEDELEC_TEST_SKIP="06,09" julia --project=. test/runtests.jl
# Cada bloque imprime QUÉ comprueba, POR QUÉ el valor esperado
# es el correcto y las magnitudes medidas, además del pasa/falla.
const TEST_VERBOSE = lowercase(get(ENV, "NEDELEC_TEST_VERBOSE", "")) in ("1", "true", "yes")
explica(args...) = (TEST_VERBOSE && println("  ▸ ", args..., "\n"); nothing)

# Lista de tests a omitir vía NEDELEC_TEST_SKIP. Acepta números (6, 06) o
# nombres (test_06), separados por comas o espacios; vacío = se ejecutan todos.
const TEST_SKIP = split(get(ENV, "NEDELEC_TEST_SKIP", ""), r"[,\s]+"; keepempty=false)
function skip_test(f)
    for s in TEST_SKIP
        n = tryparse(Int, s)
        pat = n === nothing ? s : "test_" * lpad(n, 2, '0') * "_"
        occursin(pat, f) && return true
    end
    return false
end

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
# 08     ensamblaje global (1 → 2 → 8 elem) + autovalores nulos (núcleo de K)
#        + comparación de K/M globales MATLAB vs Julia (se salta si faltan CSV)
# 09     etiquetas de la malla: contorno (PolyFace) y material (Solid)
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
    "test_09_tags.jl",
    "test_10_half_filled_cavity.jl",
)
    skip_test(f) && continue
    TEST_VERBOSE && println("\n", "═"^70, "\n▶ Ejecutando ", f, "\n", "═"^70)
    include(f)
end
