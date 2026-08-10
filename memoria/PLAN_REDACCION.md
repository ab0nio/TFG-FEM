# Plan de redacción de la memoria — entrega 4 de septiembre de 2026

Hoy: **domingo 2 de agosto de 2026**. Quedan **33 días**.

Premisa de partida: **el trabajo técnico está hecho**. El paquete Julia funciona, la suite
de 10 ficheros de test está en verde, la verificación K/M contra MATLAB está cerrada, el
refactor de etiquetas está terminado y los benchmarks están medidos. Lo que queda es
**redactar, generar figuras y escribir los dos capítulos obligatorios de la EPS que están
a cero**. Este plan asume ~5 h/día. Si algún día no llegas, recorta de las figuras
"bonitas", nunca del capítulo de verificación.

---

## Riesgo número uno: agosto

El tutor puede estar de vacaciones. **Esto se resuelve hoy, no el 25 de agosto.**
Manda un correo hoy mismo preguntando:

1. Disponibilidad en agosto y hasta qué fecha puede leer un borrador completo.
2. Si el 4 de septiembre es depósito de memoria o también defensa.
3. Extensión esperada (la EPS suele moverse en 50–70 páginas de cuerpo).
4. Si quiere memoria en español o en inglés.
5. Si acepta que el capítulo de verificación sea el núcleo del trabajo (lo es).

Hasta que responda, sigue el plan tal cual: no hay ninguna tarea bloqueada por él salvo
la revisión del borrador de la semana 4.

---

## Orden de escritura (y por qué ese orden)

**No escribas la introducción primero.** Se escribe cuando ya sabes qué has hecho. El
orden es de lo que más dominas a lo que menos:

1. **Métodos** — lo tienes en el código y en CLAUDE.md; sale casi solo.
2. **Resultados** — los números ya están medidos; es transcribir y explicar.
3. **Estado del arte** — requiere leer; por eso va después, con el trabajo ya contado.
4. **Introducción** — se escribe sabiendo ya a dónde llega el trabajo.
5. **Marco regulador y entorno socio-económico** — formularios, pero hay que documentarse.
6. **Conclusiones** — responden una a una a las preguntas de la introducción.
7. **Resumen, título y palabras clave** — **lo último**, como dice la guía de la Biblioteca.

---

## Calendario

### Domingo 2 de agosto (hoy) — 2 h

- [ ] Correo al tutor (ver arriba).
- [ ] Rellenar la portada de `memoria_tfg.tex`: titulación, curso 2025-2026, tutor, Leganés.
- [ ] Decidir el índice definitivo de capítulos y secciones, y dejarlo escrito como
      esqueleto de `\section{}` vacías. Es la tarea más rentable del mes: a partir de aquí
      escribes rellenando huecos, no mirando una página en blanco.
- [ ] Crear `memoria/figuras/` y `memoria/notas/`.

---

### Semana 1 (3–9 agosto) — MÉTODOS + los datos que ya tienes

**Objetivo: 20 páginas escritas y 6 figuras generadas.**

| Día | Tarea | Salida |
|---|---|---|
| Lun 3 | Cap. Métodos §1: formulación. Ecuación de ondas vectorial, forma débil, problema de autovalores generalizado. Leer antes: Jin cap. de elementos vectoriales + el blog de COMSOL sobre forma débil | 4 pág. |
| Mar 4 | Cap. Métodos §2: el elemento. Espacio H₂, restricción de Nédélec, los 54 monomios, tabla de DOFs (aristas/caras/volumen), G·W=I | 5 pág. |
| Mié 5 | **Figuras conceptuales** (Manim/TikZ): hexaedro de referencia con numeración de nodos, aristas y caras; distribución de los 54 DOFs; heatmaps de G y W | 3 fig. |
| Jue 6 | Cap. Métodos §3: elemento real. Mapeo isoparamétrico, jacobiano, **Piola covariante y contravariante**. Leer antes: Monk cap. correspondiente | 4 pág. |
| Vie 7 | Cap. Métodos §4: ensamblaje global. Numeración de DOFs, orientación canónica mín→máx de aristas, caras por clave ordenada, DOFs de volumen. Leer antes: Agelek 2018 (corto) | 4 pág. |
| Sáb 8 | Cap. Métodos §5: **etiquetas de malla**. PolyFace/Solid, por qué etiqueta y no geometría. Este apartado es tuyo y es bueno: el argumento del bug de εᵣ va aquí | 3 pág. |
| Dom 9 | Colchón / repaso de lo escrito | — |

**Figuras de esta semana** (datos en `visualization/manim/data/` y `NedelecHex/output/`):

- `fig_hexaedro_referencia` — nodos, aristas, caras numerados según `ReferenceHex.NODES`
- `fig_54_dofs` — escena Manim `CavidadSemillena54DOFs` ya existente, exportar frame
- `fig_matriz_G` y `fig_matriz_W` — desde `G_matrix.csv` / `W_matrix.csv`
- `fig_continuidad_tangencial` — desde `face_continuity.csv` (máx |N_A−N_B| ≈ 8e-14)
- `fig_pipeline` — diagrama de flujo malla → G/W → ensamblaje → Dirichlet → eigs

---

### Semana 2 (10–16 agosto) — RESULTADOS + ESTADO DEL ARTE

**Objetivo: 18 páginas más y las 6 figuras de datos.**

| Día | Tarea | Salida |
|---|---|---|
| Lun 10 | Cap. Resultados §1: verificación del elemento. Tabla II del paper reproducida, dualidad gᵢ(Nⱼ)=δᵢⱼ, números de condición 19 y 35 | 3 pág. |
| Mar 11 | Cap. Resultados §2: **el hallazgo K/M**. La tabla de difrel por malla, K bit-idéntica, M no, aislar el material, cerrar con el εᵣ de MATLAB. Es tu mejor resultado: dedícale espacio | 5 pág. |
| Mié 12 | **Figuras de datos**: convergencia O(h⁴), espectro con el clúster espurio λ≈0, malla coloreada por material con corte en x vs. z (la figura del bug) | 3 fig. |
| Jue 13 | Cap. Resultados §3: cavidad semillena. Autovalores analíticos vs. FEM, errores medios de las 5 mallas, filtrado de modos espurios | 4 pág. |
| Vie 14 | Cap. Resultados §4: rendimiento. Tabla por etapas Julia/MATLAB, memoria del bloque libre, la anomalía de la malla 3 declarada | 4 pág. |
| Sáb 15 | **Figuras de rendimiento**: barras logarítmicas por etapa, nnz y memoria | 3 fig. |
| Dom 16 | Colchón | — |

**Lecturas de la semana** (30–45 min cada una, tomando notas para el estado del arte):
Nédélec 1980 y 1986 (solo las definiciones), García-Castillo 2000, Amor-Martin 2016,
Amor-Martin 2019 (esta entera: es tu ensamblaje), Webb 1999, Graglia-Peterson 2011,
Davidson 2003, Falk-Gatto-Monk 2011.

**Aviso sobre honestidad en los resultados:** la varianza de los benchmarks es del ±35 %
entre tiradas. Cita los ratios redondeados (≈13×, ≈1×) y di explícitamente que son de una
máquina concreta. La anomalía de la malla 3 (0,02×) **se declara**, no se esconde: es una
decisión de corrección del umbral denso/disperso, y explicarla da más credibilidad que
ocultarla.

---

### Semana 3 (17–23 agosto) — LO QUE FALTA POR COMPLETO

**Objetivo: borrador completo en manos del tutor el domingo 23.**

| Día | Tarea | Salida |
|---|---|---|
| Lun 17 | Cap. Estado del arte: jerárquicas vs. interpolatorias/sistemáticas, alto orden en hexaedros, hueco que cubre este trabajo | 6 pág. |
| Mar 18 | Cap. Introducción: motivación, objetivos (numerados, para poder responderlos uno a uno en conclusiones), estructura del documento | 4 pág. |
| Mié 19 | **Marco regulador**. Ver guion abajo | 4 pág. |
| Jue 20 | **Entorno socio-económico**. Ver guion abajo | 4 pág. |
| Vie 21 | Cap. Conclusiones: una subsección por objetivo de la introducción, limitaciones honestas, líneas futuras | 4 pág. |
| Sáb 22 | Lista de abreviaturas, anexos, **Declaración de uso de IA generativa** (ver abajo), repaso de citas: que toda afirmación fuerte tenga referencia | — |
| **Dom 23** | **Compilar el PDF completo y enviarlo al tutor** | **HITO** |

#### Guion del marco regulador (no transcribas normativa: analiza)

- **Legislación**: el trabajo es software de simulación, no trata datos personales, así que
  el RGPD no aplica — **dilo y justifica por qué**, que es más valioso que copiar el
  articulado. Sí aplica: responsabilidad profesional del ingeniero al firmar resultados de
  simulación, y el deber de verificar el código antes de usarlo en diseño (aquí enlazas con
  tu propio capítulo de verificación y con Roache).
- **Estándares técnicos**: IEEE 754 (aritmética en punto flotante, que es lo que hace que
  ‖K‖_F coincida a 1e-15 y no exacto), el estándar del formato de malla de Gmsh, y las
  licencias de Julia (MIT), MATLAB (propietaria) y las bibliotecas usadas (ARPACK).
- **Propiedad intelectual**: el código de partida es de Amor-Martin y García-Castillo,
  publicado en Zenodo/GitHub con licencia concreta — **mírala y respétala**. Tu port a
  Julia es obra derivada: explica bajo qué licencia lo publicas y por qué es compatible.
  Este apartado te sale redondo porque es un caso real, no hipotético.

#### Guion del entorno socio-económico

- **Presupuesto**: base defendible = 12 ECTS × 25 h = **300 h** de trabajo del alumno.
  Conceptos: horas de ingeniero junior, horas de tutoría, amortización del portátil
  (fracción de 4 años), licencia de MATLAB (académica), software libre a coste 0 (Julia,
  Gmsh, LaTeX — y aquí ya tienes un argumento de impacto), costes indirectos ~15–20 %.
  El historial de git **no sirve** como justificante: solo hay 12 commits en 7 días
  distintos, no refleja las horas reales.
- **Impacto**: el resultado es una implementación libre y verificada de un elemento que
  antes solo existía en MATLAB propietario. Impacto económico = coste de licencias evitado
  en grupos de investigación y PYMEs de simulación EM. Impacto social = reproducibilidad y
  acceso. Impacto medioambiental = el ensamblaje 13× más rápido es 13× menos energía por
  simulación; si te atreves, estima kWh con el TDP del portátil, aunque sea orden de magnitud.

#### Declaración de uso de IA generativa

Has usado Claude Code de forma intensiva en este proyecto. **Decláralo con detalle y sin
adornos**: herramienta y versión, para qué exactamente (refactor de etiquetas, redacción de
tests, análisis de la discrepancia de εᵣ, generación de escenas Manim, revisión de estilo
de la memoria), qué partes del trabajo toca y **cómo verificaste lo generado** — y aquí
tienes la mejor respuesta posible: la suite de tests, la comparación numérica contra MATLAB
y los invariantes de K/M. Consúltalo con el tutor esta semana, no el 3 de septiembre.

---

### Semana 4 (24–30 agosto) — REVISIÓN

| Día | Tarea |
|---|---|
| Lun 24–Mar 25 | Mientras el tutor lee: pasada de coherencia propia. Que cada objetivo de la introducción tenga su conclusión; que cada figura esté citada en el texto y aporte algo; que cada tabla tenga fuente |
| Mié 26 | Incorporar el feedback del tutor (reserva el día entero: siempre es más de lo que esperas) |
| Jue 27 | **Resumen** (máx. 1 página: objetivo → métodos → resultados → conclusiones, sin citas, sin fórmulas, sin referencias a figuras) |
| Vie 28 | **Título definitivo** y **5 palabras clave** de la IEEE Taxonomy. Candidatas: *finite element methods*, *computational electromagnetics*, *eigenvalues and eigenfunctions*, *software verification and validation*, *scientific computing*. Rellenar `memoria_tfg.xmpdata` con los mismos datos |
| Sáb 29 | Dedicatoria/agradecimientos. Repasar que la portada no lleva NIA, DNI ni correo |
| Dom 30 | Segunda pasada de feedback si el tutor ha contestado |

---

### Semana 5 (31 agosto – 4 septiembre) — CIERRE

| Día | Tarea |
|---|---|
| Lun 31 | Revisión de formato: numeración romana I–VI correcta, páginas de cortesía en su sitio, índices regenerados, pies de figura y tabla con el estilo de ingeniería |
| Mar 1 | Ortografía y gramática de principio a fin. Comprobar que el PDF sale como **PDF/A-1b** (lo hace `pdfx`, pero verifícalo) |
| Mié 2 | Turnitin / detección de plagio si tu titulación lo exige antes del depósito. Revisar que las citas literales van entrecomilladas y referenciadas |
| Jue 3 | **Colchón intocable.** No planifiques nada. Aquí es donde aparece el problema que no viste |
| **Vie 4** | **Entrega** |

---

## Objetivos de páginas por capítulo

| Capítulo | Páginas | Estado hoy |
|---|---|---|
| 1. Introducción | 4–5 | vacío |
| 2. Estado del arte | 6–8 | vacío |
| 3. Métodos | 18–20 | **el material existe en el código** |
| 4. Resultados | 14–16 | **los números están medidos** |
| 5. Discusión y conclusiones | 4–5 | vacío |
| 6. Marco regulador | 3–4 | vacío |
| 7. Entorno socio-económico | 3–4 | vacío |
| **Cuerpo** | **52–62** | |

---

## Inventario de figuras (13)

Todas tienen ya los datos generados. Ninguna requiere volver a ejecutar el pipeline salvo
que quieras rehacer una malla.

| # | Figura | Fuente de datos |
|---|---|---|
| 1 | Hexaedro de referencia: nodos, aristas, caras | `ReferenceHex.jl` |
| 2 | Los 54 DOFs en su posición | escena Manim `CavidadSemillena54DOFs` |
| 3 | Matriz G (heatmap) | `visualization/manim/data/G_matrix.csv` |
| 4 | Matriz W = G⁻¹ (heatmap, coeficientes enteros) | `W_matrix.csv` |
| 5 | Continuidad tangencial entre dos elementos | `face_continuity.csv` |
| 6 | Diagrama de flujo del pipeline completo | — (TikZ) |
| 7 | Geometría de la cavidad semillena + malla | `matlab_nodes_*.csv` |
| 8 | Permutación de nodos de Gmsh: jacobiano antes/después | `test_06` |
| 9 | Convergencia O(h⁴) de los autovalores | `convergence_study_matlab` |
| 10 | Espectro completo: clúster espurio λ≈0 vs. modos físicos | `output/eig_all_mesh1_julia.txt` |
| 11 | **Malla coloreada por εᵣ: corte en x vs. corte en z** | `matlab_eps_*.csv` |
| 12 | Tiempos por etapa Julia vs. MATLAB (barras log) | `output/results_{julia,matlab}.csv` |
| 13 | nnz y memoria del bloque libre-libre | `output/benchmark_report.md` |

La 11 es la más valiosa del trabajo: enseña visualmente un bug real que solo se detectó
comparando matrices. Dedícale una figura grande y bien explicada.

---

## Plan de lectura (mínimo viable)

No intentes leer las 42 referencias. Estas son las que hay que leer de verdad, repartidas
en las semanas 1 y 2:

**Leer entero:** Amor-Martin 2023 (el paper, otra vez y con lápiz), Amor-Martin 2019
(ensamblaje), Agelek 2018 (orientación de aristas).

**Leer los capítulos que te tocan:** Jin 2015 (elementos vectoriales + autovalores),
Monk 2003 (H(curl) y Piola), Reddy 1993 (forma débil, si necesitas repasar la base).

**Leer el resumen y las conclusiones, para citar con criterio:** Nédélec 1980/1986,
García-Castillo 2000, Amor-Martin 2016, Webb 1999, Graglia-Peterson 2011, Davidson 2003,
Falk 2011, Boffi 2010, Roache 1998, García-Doñoro 2016.

**Consultar solo cuando escribas ese párrafo concreto:** el resto.

---

## Si vas mal de tiempo (versión mínima)

Recorta en este orden, de lo primero que se corta a lo último:

1. Figuras 6, 8 y 13 (prescindibles).
2. Estado del arte de 8 a 4 páginas.
3. Rendimiento de 4 a 2 páginas (una tabla y un párrafo).

**Nunca recortes:** el capítulo de verificación (§ Resultados 2), el marco regulador ni el
entorno socio-económico. Los dos últimos son obligatorios en la matriz de evaluación de la
EPS y su ausencia se penaliza directamente; el primero es lo que hace que tu TFG sea un
trabajo de ingeniería y no un ejercicio de programación.
