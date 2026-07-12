using NedelecHex
using Test
using LinearAlgebra, SparseArrays, Printf
# Para ejecutar:  julia --project=. test/runtests.jl
# Modo verboso:   NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl
#                 Cada bloque imprime QUÉ comprueba, POR QUÉ el valor esperado
#                 es el correcto y las magnitudes medidas, además del pasa/falla.
const TEST_VERBOSE = lowercase(get(ENV, "NEDELEC_TEST_VERBOSE", "")) in ("1", "true", "yes")
explica(args...) = (TEST_VERBOSE && println("  ▸ ", args...); nothing)
#
# Orden lógico: cada nivel depende solo de los anteriores.
#   01-02  infraestructura matemática (cuadratura, espacio H₂)
#   03-04  funcionales DOF y coeficientes (Tabla II del paper)
#   05     elemento de referencia (Jacobiano, base dual, Kᵉ/Mᵉ)
#   06     malla (permutación de nodos de Gmsh)
#   07-08  ensamblaje global (08 = dos elementos adyacentes, pedido explícito)
#   09-10  verificación física (MMS y cavidad semillena)

include("test_01_quadrature.jl")
include("test_02_monomials.jl")
include("test_03_dofs.jl")
include("test_04_coefficients.jl")
include("test_05_reference_element.jl")
include("test_06_gmsh_node_order.jl")
include("test_07_assembly.jl")
include("test_08_two_elements.jl")
include("test_09_mms.jl")
include("test_10_half_filled_cavity.jl")
