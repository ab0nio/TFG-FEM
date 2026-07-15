# TFG: Cálculo de Funciones de Base Hexahédricas de Elementos Finitos para Electromagnetismo en Julia

Trabajo de Fin de Grado (Ingeniería Telecomunicación, UC3M) sobre el Método de
Elementos Finitos aplicado a electromagnetismo computacional: implementación enn Julia
de funciones de base hexaédricas *curl-conforming* de Nédélec de segundo orden descrito en:

> A. Amor-Martin y L. E. Garcia-Castillo, "Second-Order Nédélec Curl-Conforming
> Hexahedral Element for Computational Electromagnetics," *IEEE Transactions on
> Antennas and Propagation*, vol. 71, no. 1, enero 2023.
> DOI: [10.1109/TAP.2022.3216554](https://doi.org/10.1109/TAP.2022.3216554)

## Estructura del repositorio

- **[`NedelecHex/`](NedelecHex/)**: paquete Julia con la implementación con la siguiente estructura en el orden que recomiendo leerlos para entender el proyecto:

Geometry/ReferenceHex.jl: Define la Topología del Hexahedro de referencia unidad.
Math/GaussQuadrature.jl: Calcula la integral mediante Gauss para una arista cara y volumen 
Geometry/Mesh.jl: Construye una malla estructurada y orientada a base de hexahedros(se usa sobre todo el struct que define , ya que la malla la importamos del paper).
Bases/Monomials.jl: Genera los 54 monomios con las restricciones de Nedelec y tambien el curl de los monomios.
Bases/DOFs.jl: Implementa los 3 tipos de funcionales para un monomio dado y devuelve un número real.
Bases/CoefficientGenerator.jl: Construye G(monomios evaluados) y la invierte para sacar W(coeficientes) y comprueba que tenga un número de condicionamiento bajo
Math/CoordinateTransforms.jl: Transforma un punto 2D de una cara y devuelve el correspondiente 3D
Bases/BasisEval.jl: Evalúa las 54 funciones de Base en un Punto
Assembly/GlobalDofs.jl: Pasa de DOFs locales a DOFs Globales(que dos elemntos que compartan arista tengan en esa arista mismo DOF global)
Assembly/Assembly.jl: Toma las funciones de base locales , la tabla de conectividad de los DOFs Globales y produce el sistema global Kx = \lambda Mx
Benchmark/HalfFilledCavity.jl: Cavidad Rectangular a*b*c con paredes PEC y una losa dieléctrica con permitividad /epsilon r que ocupa z perteneciente [0,h] y el resto relleno de /epsilon 0 es vacio. Se buscan los k para los que que existe un campo no trivial que satisface Maxwell , PEC en las paredes y continuidad de las componentes tangenciales resolviendo la ecuacion de resonancia
NedelecHex.jl: Módulo con todos los includes


- **`matlab/`** *(EL repo original del paper es: https://github.com/aamorm/SystematicBasesForHexahedron?tab=readme-ov-file)* 
Usado para comparar resultados.

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
| `test_01_quadrature.jl` |Comprobar la integración numeŕica por Gauss en  1D/2D/3D |
| `test_02_monomials.jl` | LComprobar que los 54 monomios y sus rotacionales están bien generados |
| `test_03_dofs.jl` | Comprobar que los integrales de los funcionales de los DOFs de Arista, Cara y Volumen se calculan correctamente comparando con integrales conocidas |
| `test_04_coefficients.jl` |Comprobar que la matriz G es unisolvente, W con coeficientes enteros y menores a 216  y coinciden con la Tabla II del paper |
| `test_05_reference_element.jl` | Jacobiano, propiedad dual, Kᵉ/Mᵉ, Tabla III (números de condición) |
| `test_06_gmsh_node_order.jl` | Comprobar que un cambio que una permutación que apliqué debido a un comportamiento inesperado de Gmsh es realmente necesario|
| `test_07_bases_eval.jl` | Las 54 bases evaluadas en el punto central, comparadas con la Tabla II del paper y con la salida de MATLAB (equivalente a `mainEvaluateCoefficients.m`) |
| `test_08_assembly.jl` | Ensamblaje global en progresión 1 → 2 → 8 elementos: un elemento (núcleo espurio de dimensión 26), dos elementos adyacentes (DOFs compartidos y orientación de aristas) y malla 2×2×2 (conteo combinatorio y propiedades globales) |
| `test_09_mms.jl` | Método de Soluciones Manufacturadas: recuperación exacta y convergencia O(h²) |
| `test_10_half_filled_cavity.jl` | Benchmark de la cavidad semillena: convergencia O(h⁴) en autovalores y equivalencia 1:1 con `mainHalfFilled.m` sobre las mallas exportadas de MATLAB |

## Verificación

- **Benchmark de la cavidad semillena**: problema de autovalores generalizado
  con solución analítica conocida (ecuación de resonancia trascendente);
  convergencia O(h²) esperada para el elemento de segundo orden.
- **Método de Soluciones Manufacturadas (MMS)**: recuperación exacta de una
  solución dentro de H₂ y convergencia log-log con pendiente ≈ 2 para una
  onda plana oblicua.
- **Cruce con MATLAB**: comparación  de matrices locales y coeficientes
  entre la implementación Julia y la referencia MATLAB.
