# TFG: Elemento Hexaédrico de Nédélec de Segundo Orden

Trabajo de Fin de Grado (Ingeniería Telecomunicación, UC3M) sobre el Método de
Elementos Finitos aplicado a electromagnetismo computacional: implementación
del elemento hexaédrico *curl-conforming* de Nédélec de segundo orden descrito en

> A. Amor-Martin y L. E. Garcia-Castillo, "Second-Order Nédélec Curl-Conforming
> Hexahedral Element for Computational Electromagnetics," *IEEE Transactions on
> Antennas and Propagation*, vol. 71, no. 1, enero 2023.
> DOI: [10.1109/TAP.2022.3216554](https://doi.org/10.1109/TAP.2022.3216554)

## Estructura del repositorio

- **[`NedelecHex/`](NedelecHex/)** — paquete Julia con la implementación
  activa: geometría del elemento de referencia, cuadratura de Gauss, los 54
  monomios y funcionales DOF del espacio H₂, generación de coeficientes,
  evaluación de bases, ensamblaje global y el benchmark de la cavidad
  semillena. Ver [`CLAUDE.md`](CLAUDE.md) para el mapeo detallado entre el
  paper y el código.
- **`matlab/`** *(no versionado, solo local)* — implementación MATLAB de
  referencia de los autores del paper, usada para verificación numérica
  cruzada frente a `NedelecHex/`.
- **[`visualization/manim/`](visualization/manim/)** — animaciones en
  Python/Manim que ilustran los 54 DOFs, la matriz de coeficientes G→W y la
  continuidad tangencial entre elementos, generadas con datos reales
  exportados de Julia.

## Uso rápido (paquete Julia)

```bash
cd NedelecHex
julia --project=.

# Suite de tests completa
julia --project=. test/runtests.jl

# Modo verboso: cada bloque explica qué comprueba y por qué
NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl
```

La suite de tests (`test/runtests.jl`) está ordenada por dependencias, desde
la cuadratura y el espacio polinómico H₂ hasta el ensamblaje global y el
benchmark físico de la cavidad semillena con eigenvalores analíticos:

| Test | Qué verifica |
|---|---|
| `test_01_quadrature.jl` | Reglas de Gauss 1D/2D/3D |
| `test_02_monomials.jl` | Los 54 monomios de H₂ y sus rotacionales analíticos |
| `test_03_dofs.jl` | Funcionales DOF de arista/cara/volumen |
| `test_04_coefficients.jl` | Matriz G unisolvente, W entero y coincide con la Tabla II del paper |
| `test_05_reference_element.jl` | Jacobiano, propiedad dual, Kᵉ/Mᵉ, Tabla III (números de condición) |
| `test_06_gmsh_node_order.jl` | Permutación de nodos de Gmsh vs. convención del elemento de referencia |
| `test_07_assembly.jl` | Ensamblaje de un elemento y de una malla 2×2×2 |
| `test_08_two_elements.jl` | DOFs compartidos entre elementos adyacentes |
| `test_09_mms.jl` | Método de Soluciones Manufacturadas: recuperación exacta y convergencia O(h²) |
| `test_10_half_filled_cavity.jl` | Benchmark de la cavidad semillena: convergencia O(h⁴) en autovalores |

## Verificación

- **Benchmark de la cavidad semillena** — problema de autovalores generalizado
  con solución analítica conocida (ecuación de resonancia trascendente);
  convergencia O(h²) esperada para el elemento de segundo orden.
- **Método de Soluciones Manufacturadas (MMS)** — recuperación exacta de una
  solución dentro de H₂ y convergencia log-log con pendiente ≈ 2 para una
  onda plana oblicua.
- **Cruce con MATLAB** — comparación 1:1 de matrices locales y coeficientes
  entre la implementación Julia y la referencia MATLAB de los autores.
