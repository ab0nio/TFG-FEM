# Ideas de visualizaciones Manim para el TFG (NedelecHex)

Criterio de coste: el portátil (Ryzen 5600H, GPU integrada, renderer Cairo) renderiza
bien escenas con **VMobjects planos** (`Line`, `Dot`, `Polygon`, `MathTex`) incluso en 3D,
pero sufre con **mallas 3D reales** (`Surface`, `Arrow3D`, muchos `Line3D`/`Dot3D`) y con
rotaciones ambientales largas (cada frame re-proyecta todo). Regla práctica:
- ✅ barato: escenas 2D, wireframes con `Line`, heatmaps como `ImageMobject`
- ⚠️ medio: `Line3D`/`Dot3D` en decenas, rotación de cámara corta
- ❌ caro: `Surface` con resolución alta, `Arrow3D` en masa, `ambient_camera_rotation` > 5 s

> **Convención de nombres**: cada escena debe exportar un PNG con el nombre de la
> sección de la memoria a la que ilustra (p. ej. `piola.png`, `convergencia.png`),
> y guardarse en `memoria/figuras/`. La tabla de abajo indica en qué sección de
> `memoria_tfg.tex` debe insertarse cada imagen (ver el comentario `% FIGURAS`
> al final de la memoria).

## Implementadas

1. **`cavidad_54dofs.py`** — pipeline completo: cavidad semillena → mallado 4×1×4 →
   zoom a un hexaedro → 54 DOFs en sus posiciones → integrales funcionales (ecs. 6–8)
   → construcción de G y W = G⁻¹ con las matrices reales de Julia. (coste medio)
   → **Memoria: [`sec:dofs`](memoria/memoria_tfg.tex:796) y [`sec:coeficientes`](memoria/memoria_tfg.tex:808)** (Métodos).
2. **`base_dual_tabla2.py`** — versión standalone y más detallada del paso G → W = G⁻¹,
   con heatmaps de las matrices reales exportadas por `CoefficientGenerator.jl`. (barato)
   → **Memoria: [`sec:coeficientes`](memoria/memoria_tfg.tex:808) y [`sec:res-tabla2`](memoria/memoria_tfg.tex:935)** (Métodos + Resultados).
3. **`continuidad_tangencial.py`** — dos elementos que comparten cara; los vectores de la
   función base evaluados desde ambos lados coinciden (datos reales de `BasisEval.jl`,
   max|N_A−N_B| ≈ 8e-14). (coste medio)
   → **Memoria: [`sec:hcurl`](memoria/memoria_tfg.tex:722)** (Estado del arte).
4. **`regularidad.py`** — por qué la forma fuerte exige C² y el FEM solo da C⁰:
  panel superior con la solución exacta suave frente a la aproximación FEM
  poligonal a trozos (codos en las fronteras de elemento), y panel inferior con
  sus derivadas (la FEM salta en cada frontera). 2D y muy barato.
  → **Memoria: [`subsec:mef-sistema`](memoria/memoria_tfg.tex:699)** (Métodos),
  figura [`fig:regularidad`](memoria/memoria_tfg.tex:713).

## Pendientes (ordenadas por relación valor/coste)

4. **Transformada de Piola covariante** (barato) — hexaedro de referencia → hexaedro
   físico distorsionado; un vector tangente a una arista se transforma con J⁻ᵀ y sigue
   siendo tangente. Es *el* concepto que justifica la conformidad H(curl). Datos reales:
   `piola_transform` de `BasisEval.jl` sobre el hexaedro distorsionado ya usado en tests.
   → **Memoria: [`sec:piola`](memoria/memoria_tfg.tex:828)** (Métodos). **ALTA prioridad**: es el
   concepto central del trabajo y no hay figura que lo ilustre.
5. **Convergencia O(h²)** (muy barato, 2D) — gráfica log-log animada con los datos reales
   de `convergence_study` (cavidad) y `test_mms_convergence.jl` (MMS, pendiente ≈ 2.0).
   Cada punto aparece con su malla en miniatura al lado.
   → **Memoria: [`sec:res-cavidad`](memoria/memoria_tfg.tex:1000)** (Resultados). **ALTA prioridad**:
   es la figura que cierra el capítulo de resultados.
6. **La permutación de Gmsh [3,7,8,4,2,6,5,1]** (barato) — numeración transfinita de Gmsh
   vs convenio `ReferenceHex.NODES`; sin permutar, el Jacobiano sale rotación pura en vez
   de diagonal (mostrar ambas matrices). Historia real de bug del proyecto → buen material
   de defensa.
   → **Memoria: [`sec:orden-nodos`](memoria/memoria_tfg.tex:843)** (Métodos).
7. **Pipeline MMS** (muy barato, 2D) — diagrama animado: u_MMS manufacturada → f = ∇×∇×u − k²u
   → resolver con FEM → comparar → pendiente 2. Con los números reales de `test_mms.jl`
   (error ~1e-13 en exactitud H₂).
   → **Memoria: [`sec:metodologia-verificacion`](memoria/memoria_tfg.tex:915)** (Métodos).
8. **Modos espurios / núcleo de K** (barato, 2D) — por qué K tiene núcleo enorme
   (gradientes ⇒ ∇×∇φ = 0), qué le hace eso a shift-invert de ARPACK, y por qué el
   benchmark usa el solver denso. Espectro real: cluster de ~0 + modos físicos.
   → **Memoria: [`sec:nucleo`](memoria/memoria_tfg.tex:896) y [`sec:espurios`](memoria/memoria_tfg.tex:737)**.
   **ALTA prioridad**: conecta la teoría de modos espurios con la resolución numérica.
9. **Campo del modo resonante** (⚠️ medio) — flechas 2D del autovector FEM en el plano
   medio y=b/2 de la cavidad, comparando k_c FEM vs analítico. Usar `Line`+`Dot` planos,
   no `Arrow3D`.
   → **Memoria: [`sec:res-cavidad`](memoria/memoria_tfg.tex:1000)** (Resultados).
10. **Cuadratura de Gauss en el hexaedro** (barato) — los 27 puntos (regla 3³) como
    `Dot3D` con radio ∝ peso; mostrar que integra exactamente polinomios de grado ≤ 5.
    → **Memoria: [`sec:KM`](memoria/memoria_tfg.tex:866)** (Métodos).
11. **Ensamblaje y patrón de dispersión** (barato, 2D) — spy-plot de K creciendo elemento
    a elemento (heatmap `ImageMobject` actualizado); DOFs de arista/cara compartidos se
    suman entre elementos.
    → **Memoria: [`sec:ensamblaje`](memoria/memoria_tfg.tex:858)** (Métodos).
12. **Orientación de aristas globales** (barato) — dos elementos vecinos cuya arista local
    va en sentidos opuestos; sin el swap de `GlobalDOFs.jl` los pesos (1−t)/t quedan
    cruzados. Otro bug-real → material de defensa.
    → **Memoria: [`sec:dofs-globales`](memoria/memoria_tfg.tex:860)** (Métodos).
13. **Estructura tensorial de H₂ / por qué 54** (muy barato, 2D) — tabla animada de grados
    por componente (ξ: 1×2×2, …), restricción de Nédélec de orden mixto, 3×18 = 54.
    → **Memoria: [`sec:h2`](memoria/memoria_tfg.tex:788)** (Métodos).

## Resumen de prioridad

| Prioridad | Escena | Sección memoria |
|---|---|---|
| 🔴 Alta | 4. Piola covariante | [`sec:piola`](memoria/memoria_tfg.tex:828) |
| 🔴 Alta | 5. Convergencia O(h²) | [`sec:res-cavidad`](memoria/memoria_tfg.tex:1000) |
| 🔴 Alta | 8. Modos espurios / núcleo de K | [`sec:nucleo`](memoria/memoria_tfg.tex:896) |
| 🟡 Media | 1. Cavidad 54 DOFs | [`sec:dofs`](memoria/memoria_tfg.tex:796) |
| 🟡 Media | 2. Base dual / Tabla II | [`sec:res-tabla2`](memoria/memoria_tfg.tex:935) |
| 🟡 Media | 3. Continuidad tangencial | [`sec:hcurl`](memoria/memoria_tfg.tex:722) |
| 🟡 Media | 9. Campo del modo resonante | [`sec:res-cavidad`](memoria/memoria_tfg.tex:1000) |
| 🟢 Baja | 6. Permutación Gmsh | [`sec:orden-nodos`](memoria/memoria_tfg.tex:843) |
| 🟢 Baja | 7. Pipeline MMS | [`sec:metodologia-verificacion`](memoria/memoria_tfg.tex:915) |
| 🟢 Baja | 10. Cuadratura Gauss | [`sec:KM`](memoria/memoria_tfg.tex:866) |
| 🟢 Baja | 11. Ensamblaje / spy-plot | [`sec:ensamblaje`](memoria/memoria_tfg.tex:858) |
| 🟢 Baja | 12. Orientación de aristas | [`sec:dofs-globales`](memoria/memoria_tfg.tex:860) |
| 🟢 Baja | 13. Estructura tensorial H₂ | [`sec:h2`](memoria/memoria_tfg.tex:788) |

> **Nota**: las escenas 1–3 figuran como "implementadas" en el listado original, pero
> los ficheros `base_dual_tabla2.py` y `continuidad_tangencial.py` no están en el
> directorio. Verificar si existen o están pendientes de crear.
