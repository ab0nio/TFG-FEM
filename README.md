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


- **`matlab/`** *(EL repo original del paper es: https://github.com/aamorm/SystematicBasesForHexahedron?tab=readme-ov-file)*
Implementación de referencia del paper, usada para comparar resultados. Usa clases OOP (`@Clase`) y espacios de nombres (`+paquete`).:

  **Construir los monomios**
  1. **`+math/@RealCoordinates/RealCoordinates.m`**: Es una clase que guarda las variables simbólicas x, y, z (no valores concretos, sino literalmente los símbolos (x, y ,z) para crear funciones), y que se pasa por referencia porque hereda de handle.
  2. **`+math/@Polynomial/Polynomial.m`**: Pasa de una funcón simbólica a un array de coeficientes para operar más rápido
  3. **`+math/@Polynomial/evaluatePolynomial.m`**: evalúa un polinomio en un punto r.
  4. **`+math/@Polynomial/getSymPolynomial.m`**: devuelve la forma simbólica del polinomio .
  5. **`+math/getMixedPolySpace.m`**: para un array t de variables simbólicas y un determinado orden para cada variable construye el espacio polinómico `Q(o1,o2,o3)` haciendo el producto tensorial de las 3 1D.
  6. **`+math/getFullPolySpace.m`**: devuelve el vector {1, t, t², …, tⁿ} con todas las potencias de una única variable hasta el grado n.
  7. **`+bases/getPolynomials.m`**: genera los **54 monomios** como producto tensorial de los espacios mixtos (equivalente a `Monomials.jl`).

  **Cuadratura de Gauss**
  8. **`+math/getGaussPointsForSegment.m`**: pesos y puntos de Gauss 1D (para aristas).
  9. **`+math/getGaussPointsForQuadrilateral.m`**: Gauss 2D en una cara (producto tensorial del 1D).
  10. **`+math/getGaussPoints.m`**: Gauss 3D en el hexaedro (producto tensorial 2D*1D). 
  Equivalen a `GaussQuadrature.jl`.
  11. **`+math/transform2Dto3DInReferenceElement.m`**: pasa un punto 2D de una cara a su 3D en el elemento de referencia. Equivale a `CoordinateTransforms.jl`.

  **Construcción de la geometría**
  12. **`+geometry/@Node/Node.m`**: un nodo(vértice) almacena: coordenadas, identificador global y si tiene o no condición de contorno.
  13. **`+geometry/@Edge/Edge.m`** : una arista son 2 nodos, sus DOFs y CC.
  14. **`+geometry/@Face/Face.m`**: una cara sus aristas, sus DOFs y su orden 
  15. **`+geometry/@Element/Element.m`** un elemento tiene nodos+aristas+caras+DOFs interiores de volumen,material ,id ,el orden yfabrica el elemento de referencia 
  16. **`+geometry/@PolyFace/PolyFace.m`**: grupos de caras que definen las condiciones de contorno en este caso dirichlet.
  17. **`+geometry/@Solid/Solid.m`**: grupos de elementos que comparten material.
  18. **`+geometry/@Mesh/Mesh.m`**: Contiene referencias a todo lo anterior + settings del problema. No construye la malla ,pero si la carga de las del .mat(equivale al `struct HexMesh`, pero en OOP).

  **Materiales**
  19. **`+settings/@MaterialProperty/`** (`MaterialProperty.m`, `evaluateProperty.m`), **`@Material/Material.m`**, **`@SystemSettings/SystemSettings.m`**: ε_r, μ_r (posiblemente inhomogéneos/anisótropos) y la configuración del problema.

  **Las funciones base (el corazón — Fig. 2 y 3 del paper)**
  20. **`+bases/@SystematicBases/SystematicBases.m`**:construye los 54 monomios, sus rotacionales y los polinomios de test `q` de la Tabla I (o los carga si ya están calculados).
  21. **`.../evaluateEdgeDOF.m`**: funcional de arista, ec. (6).
  22. **`.../evaluateFaceDOF.m`**: funcional de cara, ec. (7).
  23. **`.../evaluateVolumeDOF.m`**: funcional de volumen, ec. (8).
  24. **`.../getCoefficients.m`**: construye la matriz **G** (DOF `i` aplicado a monomio `j`) y resuelve **W = G\I**  Equivalente a `CoefficientGenerator.jl`.
  25. **`.../getBasesEvaluated.m`**: evalúa las 54 bases en el elemento de **referencia** en un punto (combinación de monomios con W). Equivale a `eval_bases_ref` de `BasisEval.jl`. El mapeo de Piola al elemento real se hace luego, en el ensamblaje.
  26. **`.../getCurlBasesEvaluated.m`**: lo mismo pero para los **rotacionales** de las 54 bases.
  27. **`.../saveCoefficients.m`**, **`.../showBases.m`**, **`+bases/showPolynomials.m`**: guardar los coeficientes para no recalcularlos y mostrar coeficientes/monomios de forma legible.

  **Numeración global de DOFs**
  28. **`@Edge/generateDOF.m`**, **`@Face/generateDOF.m`**, **`@Element/generateVolumeDOF.m`**: asignan 2 / 4 / 6 DOFs por arista / cara / elemento.
  29. **`@Mesh/generateDOF.m`**: recorre todas las entidades y da un identificador único a cada DOF (equivale a `build_dof_map` de `GlobalDOFs.jl`).
  30. **`@Element/getDOF.m`**: lista los DOFs de un elemento en orden; **`@Mesh/getNNZ.m`**: es una cota superior de no-ceros para preasignar las dispersas.

  **Ensamblaje FEM**
  31. **`+assembly/getMassMatrix.m`**: matriz de masa del elemento `∫ εᵣ N·N` (c₂, ec. 21).
  32. **`+assembly/getStiffMatrix.m`**: matriz de rigidez del elemento `∫ (1/μᵣ) (∇×N)·(∇×N)` ( c₁, ec. 20).
  33. **`+assembly/imposeLocalDirichletMatrix.m`**: impone PEC por elemento (anula filas/columnas, 1 en la diagonal).
  34. **`+assembly/getGlobalMatrices.m`**: ensambla K y M globales recorriendo los elementos. Equivale a `assemble_global_matrices` de `Assembly.jl`.

  **Benchmark de la cavidad semillena**
  35. **`+math/getKcFromHalfFilledCavity.m`**: autovalor **analítico** usando la ecuación de resonancia trascendente. Equivale a `analytic_kc` de `HalfFilledCavity.jl`.
  36. **`@Mesh/getLengthStatistics.m`**: tamaño medio/máx/mín de arista (eje x de la gráfica de convergencia).
  37. **`+math/plotHalfFilledCavity.m`**: dibuja el error en log-log (Fig. 5).

  **Drivers (los `main` — puntos de entrada)**
  38. **`mainGenerateCoefficients.m`**: solo genera las bases (Fig. 2).
  39. **`mainEvaluateCoefficients.m`**: evalúa las bases a partir de la Tabla II (Fig. 3).
  40. **`mainHalfFilled.m`**: flujo completo (bases → mallas → ensamblaje → autovalores → error).
  41. **`mainDemo.m`** y **`mainBenchmark.m`**: comparación cruzada MATLAB↔Julia (W, bases, matrices) y benchmark de `getCoefficients`. **`exportMesh.m`**: exporta las mallas a CSV para que Julia use la misma discretización.

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
