# Benchmark de rendimiento (ensamblaje + solver): Julia vs MATLAB

## Tiempo en caliente por etapa

| Etapa | Malla | n_dof | Julia (ms) | MATLAB (ms) | MATLAB/Julia |
|---|---:|---:|---:|---:|---:|
| coefficients | — | 0 | 5.2617 | 251.0346 | 47.71× |
|  |  |  |  |  |  |
| dof_numbering | 1 | 318 | 0.0097 | 0.4276 | 44.10× |
| assembly | 1 | 318 | 8.0376 | 201.9896 † | — |
| dirichlet | 1 | 318 | 0.0363 | — | — |
| **assembly+dirichlet** | 1 | 318 | **8.0739** | **201.9896** | **25.02×** |
| solver | 1 | 318 | 0.3783 | 1.5476 | 4.09× |
|  |  |  |  |  |  |
| dof_numbering | 2 | 1678 | 0.0560 | 2.1716 | 38.79× |
| assembly | 2 | 1678 | 63.9433 | 1093.0986 † | — |
| dirichlet | 2 | 1678 | 0.1965 | — | — |
| **assembly+dirichlet** | 2 | 1678 | **64.1398** | **1093.0986** | **17.04×** |
| solver | 2 | 1678 | 50.9357 | 7.8646 | 0.15× |
|  |  |  |  |  |  |
| dof_numbering | 3 | 4062 | 0.1253 | 6.2516 | 49.88× |
| assembly | 3 | 4062 | 161.6445 | 2634.1486 † | — |
| dirichlet | 3 | 4062 | 0.5811 | — | — |
| **assembly+dirichlet** | 3 | 4062 | **162.2256** | **2634.1486** | **16.24×** |
| solver | 3 | 4062 | 1216.0006 | 19.0496 | 0.02× |
|  |  |  |  |  |  |
| dof_numbering | 4 | 14916 | 0.6391 | 18.1826 | 28.45× |
| assembly | 4 | 14916 | 927.6766 | 11361.4036 † | — |
| dirichlet | 4 | 14916 | 2.4535 | — | — |
| **assembly+dirichlet** | 4 | 14916 | **930.1301** | **11361.4036** | **12.21×** |
| solver | 4 | 14916 | 250.7212 | 260.1546 | 1.04× |
|  |  |  |  |  |  |
| dof_numbering | 5 | 23124 | 1.2037 | 28.2926 | 23.50× |
| assembly | 5 | 23124 | 1313.1491 | 17136.7146 † | — |
| dirichlet | 5 | 23124 | 4.1139 | — | — |
| **assembly+dirichlet** | 5 | 23124 | **1317.2629** | **17136.7146** | **13.01×** |
| solver | 5 | 23124 | 532.0363 | 424.2356 | 0.80× |

> Ratio > 1 → Julia más rápido; < 1 → MATLAB más rápido.
> † `assembly` **no es comparable etapa a etapa**: MATLAB impone las
> condiciones Dirichlet dentro de `getGlobalMatrices`, mientras que en
> Julia son una etapa aparte (`dirichlet`). La comparación de igual a
> igual es la fila **assembly+dirichlet** = (assembly + dirichlet) de
> Julia frente a `assembly` de MATLAB; por eso el ratio se da ahí y no
> en la fila `assembly`.

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

## Memoria de K + M — matriz completa (contexto, no comparable)

| Malla | n_dof | nnz(K) Julia | nnz(K) MATLAB | Mem Julia | Mem MATLAB |
|---:|---:|---:|---:|---:|---:|
| 1 | 318 | 21900 | 2780 | 514.8 KB | 70.4 KB |
| 2 | 1678 | 128300 | 28860 | 3.1 MB | 708.9 KB |
| 3 | 4062 | 318888 | 80284 | 7.6 MB | 1.9 MB |
| 4 | 14916 | 1321528 | 693864 | 31.3 MB | 28.2 MB |
| 5 | 23124 | 2061344 | 1108104 | 48.9 MB | 44.3 MB |

> **Estas columnas NO son comparables entre sí.** MATLAB impone Dirichlet
> DENTRO de `getGlobalMatrices` (`imposeLocalDirichletMatrix` anula fila y
> columna y pone un 1 en la diagonal), así que su matriz "completa" ya
> llega con los acoplamientos de contorno a cero — `sparse` los descarta —
> y su nnz es prácticamente el del bloque libre. La K de Julia es la matriz
> íntegra n_dof x n_dof y el filtrado es posterior. MATLAB parece gastar
> menos memoria porque guarda menos matriz, no porque la guarde mejor.
> La comparación buena es la tabla anterior.

## TTFX — primera ejecución (incluye compilación)

| Etapa | Malla | Julia TTFX (ms) | MATLAB TTFX (ms) |
|---|---:|---:|---:|
| coefficients | — | 2955.3417 | 486.0640 |
|  |  |  |  |
| dof_numbering | 1 | 0.0233 | 12.4480 |
| assembly | 1 | 318.4896 | 294.0400 |
| dirichlet | 1 | 0.0853 | — |
| **assembly+dirichlet** | 1 | **318.5749** | **294.0400** |
| solver | 1 | 0.5148 | 145.9280 |
|  |  |  |  |
| dof_numbering | 2 | 0.0696 | 3.5420 |
| assembly | 2 | 60.8644 | 1096.1230 |
| dirichlet | 2 | 0.2355 | — |
| **assembly+dirichlet** | 2 | **61.0998** | **1096.1230** |
| solver | 2 | 61.5199 | 9.1920 |
|  |  |  |  |
| dof_numbering | 3 | 0.1510 | 6.5230 |
| assembly | 3 | 159.9911 | 2979.8210 |
| dirichlet | 3 | 0.5718 | — |
| **assembly+dirichlet** | 3 | **160.5630** | **2979.8210** |
| solver | 3 | 1148.9181 | 19.2340 |
|  |  |  |  |
| dof_numbering | 4 | 0.5727 | 19.0750 |
| assembly | 4 | 965.9595 | 11238.5500 |
| dirichlet | 4 | 2.7159 | — |
| **assembly+dirichlet** | 4 | **968.6754** | **11238.5500** |
| solver | 4 | 1610.2406 | 272.6770 |
|  |  |  |  |
| dof_numbering | 5 | 1.0499 | 29.6810 |
| assembly | 5 | 1419.5994 | 17331.3430 |
| dirichlet | 5 | 3.9109 | — |
| **assembly+dirichlet** | 5 | **1423.5103** | **17331.3430** |
| solver | 5 | 465.7025 | 441.6270 |

> En Julia el TTFX de la primera malla incluye la compilación JIT; en
> las mallas siguientes la función ya está compilada, de ahí la caída.
> Misma salvedad que arriba: la fila comparable con el `assembly` de
> MATLAB es **assembly+dirichlet**.
