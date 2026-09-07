# Benchmark de rendimiento (ensamblaje + solver): Julia vs MATLAB

- Fecha: 2026-09-06 21:35:44
- Máquina: `Linux 7.1.8-arch1-3 x86_64`
- Hilos (BLAS/cómputo): 1 en ambos lados (comparación de algoritmo serial).
- Solver de autovalores: **incluido** como etapa `solver`, para poder
  contrastar su coste con el de construir las matrices.
- Tiempos: **mediana en caliente** (BenchmarkTools / timeit). Ver TTFX aparte.

## Tiempo en caliente por etapa

| Etapa | Malla | n_dof | Julia (ms) | MATLAB (ms) | MATLAB/Julia |
|---|---:|---:|---:|---:|---:|
| coefficients | — | 0 | 4.7732 | 251.0346 | 52.59× |
|  |  |  |  |  |  |
| dof_numbering | 1 | 318 | 0.0086 | 0.4276 | 49.75× |
| assembly | 1 | 318 | 6.7801 | 201.9896 | 29.79× |
| solver | 1 | 318 | 0.3845 | 1.5476 | 4.02× |
|  |  |  |  |  |  |
| dof_numbering | 2 | 1678 | 0.0517 | 2.1716 | 41.97× |
| assembly | 2 | 1678 | 54.8406 | 1093.0986 | 19.93× |
| solver | 2 | 1678 | 45.6151 | 7.8646 | 0.17× |
|  |  |  |  |  |  |
| dof_numbering | 3 | 4062 | 0.1269 | 6.2516 | 49.28× |
| assembly | 3 | 4062 | 144.6518 | 2634.1486 | 18.21× |
| solver | 3 | 4062 | 983.1298 | 19.0496 | 0.02× |
|  |  |  |  |  |  |
| dof_numbering | 4 | 14916 | 0.5328 | 18.1826 | 34.13× |
| assembly | 4 | 14916 | 1150.0440 | 11361.4036 | 9.88× |
| solver | 4 | 14916 | 474.8820 | 260.1546 | 0.55× |
|  |  |  |  |  |  |
| dof_numbering | 5 | 23124 | 0.9610 | 28.2926 | 29.44× |
| assembly | 5 | 23124 | 2563.8155 | 17136.7146 | 6.68× |
| solver | 5 | 23124 | 430.9752 | 424.2356 | 0.98× |

> Ratio > 1 → Julia más rápido; < 1 → MATLAB más rápido.
> `assembly` es comparable 1:1: en ambos lenguajes la etapa incluye la
> imposición de las condiciones Dirichlet (PEC) dentro del ensamblaje.

## Memoria de K + M — bloque libre-libre (comparable)

| Malla | n_dof | nnz(K_f) Julia | nnz(K_f) MATLAB | Mem Julia | Mem MATLAB |
|---:|---:|---:|---:|---:|---:|
| 1 | 318 | 2556 | 2556 | 61.4 KB | 59.9 KB |
| 2 | 1678 | 27868 | 27868 | 664.6 KB | 662.4 KB |
| 3 | 4062 | 78012 | 78012 | 1.8 MB | 1.8 MB |
| 4 | 14916 | 688676 | 688744 | 15.8 MB | 15.8 MB |
| 5 | 23124 | 1100144 | 1100424 | 25.2 MB | 25.2 MB |

> Ambas columnas son `K(free,free)` y `M(free,free)`: las matrices **con
> Dirichlet ya impuesto en los dos lados**, mismo problema discreto y mismo
> tamaño (n_free x n_free). Es la única comparación de memoria válida.
> `nnz(K_f)` coincide **exactamente** en las mallas 1–3 y difiere en ~1e-4
> relativo en las 4–5 (68 y 280 entradas de 6.9e5 y 1.1e6): son términos de
> redondeo que cancelan a cero exacto en un orden de acumulación y no en el
> otro — hay decenas de miles de entradas por debajo de 1e-9·‖K_f‖_F. La
> igualdad numérica de fondo la certifica `test_08_assembly.jl` (‖K‖_F y
> traza coinciden a ~1e-15).

## Memoria de K + M — matriz completa (contexto)

| Malla | n_dof | nnz(K) Julia | nnz(K) MATLAB | Mem Julia | Mem MATLAB |
|---:|---:|---:|---:|---:|---:|
| 1 | 318 | 2780 | 2780 | 144.4 KB | 70.4 KB |
| 2 | 1678 | 28860 | 28860 | 1.3 MB | 708.9 KB |
| 3 | 4062 | 80284 | 80284 | 3.6 MB | 1.9 MB |
| 4 | 14916 | 693796 | 693864 | 26.8 MB | 28.2 MB |
| 5 | 23124 | 1107824 | 1108104 | 46.5 MB | 44.3 MB |

> Ambas implementaciones imponen Dirichlet DENTRO del ensamblaje (`dirichlet=true`
> en Julia; `imposeLocalDirichletMatrix` en MATLAB), así que la matriz "completa"
> ya llega con los acoplamientos de contorno a cero y su nnz es prácticamente el
> del bloque libre. La memoria es algo mayor que la de la tabla anterior porque la
> matriz es n_dof x n_dof (no n_free x n_free), pero la comparación entre lenguajes
> es directa. La tabla anterior (bloque libre-libre) sigue siendo la referencia.

## TTFX — primera ejecución (incluye compilación)

| Etapa | Malla | Julia TTFX (ms) | MATLAB TTFX (ms) |
|---|---:|---:|---:|
| coefficients | — | 2536.0688 | 486.0640 |
|  |  |  |  |
| dof_numbering | 1 | 0.0264 | 12.4480 |
| assembly | 1 | 410.2340 | 294.0400 |
| solver | 1 | 0.5246 | 145.9280 |
|  |  |  |  |
| dof_numbering | 2 | 0.0700 | 3.5420 |
| assembly | 2 | 52.9709 | 1096.1230 |
| solver | 2 | 46.2147 | 9.1920 |
|  |  |  |  |
| dof_numbering | 3 | 0.1349 | 6.5230 |
| assembly | 3 | 143.7643 | 2979.8210 |
| solver | 3 | 1004.2435 | 19.2340 |
|  |  |  |  |
| dof_numbering | 4 | 0.7046 | 19.0750 |
| assembly | 4 | 926.3353 | 11238.5500 |
| solver | 4 | 1473.5822 | 272.6770 |
|  |  |  |  |
| dof_numbering | 5 | 1.0774 | 29.6810 |
| assembly | 5 | 2937.2297 | 17331.3430 |
| solver | 5 | 679.6993 | 441.6270 |

> En Julia el TTFX de la primera malla incluye la compilación JIT; en
> las mallas siguientes la función ya está compilada, de ahí la caída.
> El `assembly` es comparable 1:1: en ambos lenguajes la etapa incluye la
> imposición de las condiciones Dirichlet (PEC) dentro del ensamblaje.
