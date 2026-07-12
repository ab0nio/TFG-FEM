# Ideas de visualizaciones Manim para el TFG (NedelecHex)

Criterio de coste: el portátil (Ryzen 5600H, GPU integrada, renderer Cairo) renderiza
bien escenas con **VMobjects planos** (`Line`, `Dot`, `Polygon`, `MathTex`) incluso en 3D,
pero sufre con **mallas 3D reales** (`Surface`, `Arrow3D`, muchos `Line3D`/`Dot3D`) y con
rotaciones ambientales largas (cada frame re-proyecta todo). Regla práctica:
- ✅ barato: escenas 2D, wireframes con `Line`, heatmaps como `ImageMobject`
- ⚠️ medio: `Line3D`/`Dot3D` en decenas, rotación de cámara corta
- ❌ caro: `Surface` con resolución alta, `Arrow3D` en masa, `ambient_camera_rotation` > 5 s

## Implementadas

1. **`cavidad_54dofs.py`** — pipeline completo: cavidad semillena → mallado 4×1×4 →
   zoom a un hexaedro → 54 DOFs en sus posiciones → integrales funcionales (ecs. 6–8)
   → construcción de G y W = G⁻¹ con las matrices reales de Julia. (coste medio)
2. **`base_dual_tabla2.py`** — versión standalone y más detallada del paso G → W = G⁻¹,
   con heatmaps de las matrices reales exportadas por `CoefficientGenerator.jl`. (barato)
3. **`continuidad_tangencial.py`** — dos elementos que comparten cara; los vectores de la
   función base evaluados desde ambos lados coinciden (datos reales de `BasisEval.jl`,
   max|N_A−N_B| ≈ 8e-14). (coste medio)

## Pendientes (ordenadas por relación valor/coste)

4. **Transformada de Piola covariante** (barato) — hexaedro de referencia → hexaedro
   físico distorsionado; un vector tangente a una arista se transforma con J⁻ᵀ y sigue
   siendo tangente. Es *el* concepto que justifica la conformidad H(curl). Datos reales:
   `piola_transform` de `BasisEval.jl` sobre el hexaedro distorsionado ya usado en tests.
5. **Convergencia O(h²)** (muy barato, 2D) — gráfica log-log animada con los datos reales
   de `convergence_study` (cavidad) y `test_mms_convergence.jl` (MMS, pendiente ≈ 2.0).
   Cada punto aparece con su malla en miniatura al lado.
6. **La permutación de Gmsh [3,7,8,4,2,6,5,1]** (barato) — numeración transfinita de Gmsh
   vs convenio `ReferenceHex.NODES`; sin permutar, el Jacobiano sale rotación pura en vez
   de diagonal (mostrar ambas matrices). Historia real de bug del proyecto → buen material
   de defensa.
7. **Pipeline MMS** (muy barato, 2D) — diagrama animado: u_MMS manufacturada → f = ∇×∇×u − k²u
   → resolver con FEM → comparar → pendiente 2. Con los números reales de `test_mms.jl`
   (error ~1e-13 en exactitud H₂).
8. **Modos espurios / núcleo de K** (barato, 2D) — por qué K tiene núcleo enorme
   (gradientes ⇒ ∇×∇φ = 0), qué le hace eso a shift-invert de ARPACK, y por qué el
   benchmark usa el solver denso. Espectro real: cluster de ~0 + modos físicos.
9. **Campo del modo resonante** (⚠️ medio) — flechas 2D del autovector FEM en el plano
   medio y=b/2 de la cavidad, comparando k_c FEM vs analítico. Usar `Line`+`Dot` planos,
   no `Arrow3D`.
10. **Cuadratura de Gauss en el hexaedro** (barato) — los 27 puntos (regla 3³) como
    `Dot3D` con radio ∝ peso; mostrar que integra exactamente polinomios de grado ≤ 5.
11. **Ensamblaje y patrón de dispersión** (barato, 2D) — spy-plot de K creciendo elemento
    a elemento (heatmap `ImageMobject` actualizado); DOFs de arista/cara compartidos se
    suman entre elementos.
12. **Orientación de aristas globales** (barato) — dos elementos vecinos cuya arista local
    va en sentidos opuestos; sin el swap de `GlobalDOFs.jl` los pesos (1−t)/t quedan
    cruzados. Otro bug-real → material de defensa.
13. **Estructura tensorial de H₂ / por qué 54** (muy barato, 2D) — tabla animada de grados
    por componente (ξ: 1×2×2, …), restricción de Nédélec de orden mixto, 3×18 = 54.
