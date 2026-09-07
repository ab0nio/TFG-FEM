# TFG: Cálculo de Funciones de Base Hexaédricas de Elementos Finitos para Electromagnetismo en Julia

Trabajo de Fin de Grado (Ingeniería Telecomunicación, UC3M) sobre el Método de
Elementos Finitos aplicado a electromagnetismo computacional: implementación en Julia
de funciones de base hexaédricas *curl-conforming* de Nédélec de segundo orden descrito en:

> A. Amor-Martin y L. E. Garcia-Castillo, "Second-Order Nédélec Curl-Conforming
> Hexahedral Element for Computational Electromagnetics," *IEEE Transactions on
> Antennas and Propagation*, vol. 71, no. 1, enero 2023.
> DOI: [10.1109/TAP.2022.3216554](https://doi.org/10.1109/TAP.2022.3216554)

## Estructura del repositorio

- **[`NedelecHex/`](NedelecHex/)**: paquete Julia con la implementación. En el orden en que recomiendo leer los módulos para entender el proyecto:

  1. **`Geometry/ReferenceHex.jl`**: define la topología del hexaedro de referencia unidad.
  2. **`Math/GaussQuadrature.jl`**: calcula la integral mediante Gauss para arista, cara y volumen.
  3. **`Geometry/Mesh.jl`**: construye una malla estructurada y orientada a base de hexaedros (se usa sobre todo el `struct` que define, ya que la malla la importamos del paper).
  4. **`Bases/Monomials.jl`**: genera los 54 monomios con las restricciones de Nédélec y también el curl de los monomios.
  5. **`Bases/DOFs.jl`**: implementa los 3 tipos de funcionales para un monomio dado y devuelve un número real.
  6. **`Bases/CoefficientGenerator.jl`**: construye G (monomios evaluados) y la invierte para sacar W (coeficientes), y comprueba que tenga un número de condicionamiento bajo.
  7. **`Math/CoordinateTransforms.jl`**: transforma un punto 2D de una cara y devuelve el 3D correspondiente.
  8. **`Bases/BasisEval.jl`**: evalúa las 54 funciones de base en un punto.
  9. **`Assembly/GlobalDOFs.jl`**: pasa de DOFs locales a DOFs globales (que dos elementos que compartan arista tengan el mismo DOF global en esa arista).
  10. **`Assembly/Assembly.jl`**: toma las funciones de base locales y la tabla de conectividad de los DOFs globales, y produce el sistema global `Kx = λMx`.
  11. **`Benchmark/HalfFilledCavity.jl`**: cavidad rectangular a×b×c con paredes PEC y una losa dieléctrica de permitividad εᵣ que ocupa z ∈ [0, h], con el resto (ε₀) vacío. Se buscan los k para los que existe un campo no trivial que satisface Maxwell, PEC en las paredes y continuidad de las componentes tangenciales, resolviendo la ecuación de resonancia.
  12. **`NedelecHex.jl`**: módulo que agrupa todos los `include`.


- **`matlab/`**: implementación de referencia del paper, usada para comparar
  resultados (no se versiona en este repo; el original es
  https://github.com/aamorm/SystematicBasesForHexahedron).

- **[`visualization/manim/`](visualization/manim/)**: animaciones en
  Python/Manim aún en proceso.

## Uso (paquete Julia)

```bash
cd NedelecHex
julia --project=.

# Ejecutar todos los tests
julia --project=. test/runtests.jl

# Ejecutar los tests con explicaciones de lo que van comprobando
NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl
```

Los tests (`test/runtests.jl`) van en orden construyendo uno encima del otro en orden de complejidad ascendente

| Test | Qué verifica |
|---|---|
| `test_01_quadrature.jl` |Comprobar la integración numérica por Gauss en  1D/2D/3D |
| `test_02_monomials.jl` | Comprobar que los 54 monomios y sus rotacionales están bien generados |
| `test_03_dofs.jl` | Comprobar que las integrales de los funcionales de los DOFs de Arista, Cara y Volumen se calculan correctamente comparando con integrales conocidas |
| `test_04_coefficients.jl` |Comprobar que la matriz G es unisolvente, W con coeficientes enteros y menores a 216  y coinciden con la Tabla II del paper |
| `test_05_reference_element.jl` | Jacobiano, propiedad dual, Kᵉ/Mᵉ, Tabla III (números de condición) |
| `test_06_gmsh_node_order.jl` | Comprobar que una permutación que apliqué debido a un comportamiento inesperado de Gmsh es realmente necesaria|
| `test_07_bases_eval.jl` | Las 54 bases evaluadas en el punto central, comparadas con la Tabla II del paper y con la salida de MATLAB (equivalente a `mainEvaluateCoefficients.m`) |
| `test_08_assembly.jl` | Ensamblaje global en progresión 1 → 2 → 8 elementos: un elemento (núcleo espurio de dimensión 26), dos elementos adyacentes (DOFs compartidos y orientación de aristas) y malla 2×2×2 (conteo combinatorio y propiedades globales) |
| `test_09_mms.jl` | Método de Soluciones Manufacturadas: en proceso |
| `test_10_half_filled_cavity.jl` | Benchmark de la cavidad semillena: convergencia O(h⁴) en autovalores y equivalencia 1:1 con `mainHalfFilled.m` sobre las mallas exportadas de MATLAB |

## Verificación

- **Benchmark de la cavidad semillena**: problema de autovalores generalizado
  con solución analítica conocida (ecuación de resonancia trascendente);
  convergencia O(h²) esperada para el elemento de segundo orden

- **Cruce con MATLAB**: comparación  de matrices locales y coeficientes
  entre la implementación Julia y la referencia MATLAB.
