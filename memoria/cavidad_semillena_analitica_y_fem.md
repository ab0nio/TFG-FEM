# La cavidad semillena: solución analítica y solución por FEM

Documento de estudio. Cada ecuación se **presenta y se definen todos sus símbolos antes de
usarla**. Primero la **solución analítica** (la referencia exacta) y después la **solución por
FEM** (la aproximación numérica que hay que validar).

> Las fórmulas están en LaTeX. Se ven correctamente en la vista previa de Markdown de VS Code
> (o en GitHub). Si abres este archivo en un terminal plano, las verás como texto crudo.

---

## 0. El problema físico (común a las dos soluciones)

Una cavidad rectangular metálica (paredes de conductor eléctrico perfecto, PEC) rellena a
medias por una losa dieléctrica.

```
        z = c  +--------------+   <- pared PEC (arriba)
               |    VACÍO      |      región 2:  eps = 1,      grosor (c - h)
        z = h  +--------------+   <- interfaz dieléctrico / vacío
               |  DIELÉCTRICO  |      región 1:  eps = eps_r,  grosor h
        z = 0  +--------------+   <- pared PEC (abajo)

   Dimensiones: a (en x),  b (en y),  c (en z).
```

**Nomenclatura** (símbolos usados en todo el documento):

| Símbolo | Significado | Unidad |
|---|---|---|
| $a,\,b,\,c$ | dimensiones de la cavidad en $x,y,z$ | m |
| $h$ | grosor de la losa dieléctrica ($0\le z\le h$) | m |
| $\varepsilon_r$ | permitividad relativa del dieléctrico (aquí $\varepsilon_r=2$) | — |
| $\varepsilon_0,\,\mu_0$ | permitividad y permeabilidad del vacío | — |
| $\omega$ | frecuencia angular | rad/s |
| $\mathbf{E}$ | campo eléctrico (vector) | V/m |
| $k\equiv k_c$ | número de onda en el vacío, $k=\omega/c_0$ **(lo que buscamos)** | rad/m |
| $c_0$ | velocidad de la luz, $c_0=1/\sqrt{\mu_0\varepsilon_0}$ | m/s |
| $m,\,n$ | índices enteros del modo transversal | — |

Buscamos los **números de onda de resonancia** $k_c$: los valores discretos para los que la
cavidad resuena. De $k_c$ sale la frecuencia con $f = c_0\,k_c/(2\pi)$.

---

# PARTE I — Solución analítica

Idea: resolver Maxwell exactamente aprovechando que la geometría son cajas rectangulares. Se
llega a **una** ecuación cuyas raíces son los $k_c$.

## I.1 — Ecuación de Helmholtz (punto de partida)

**Se expone.** En cada región, sin fuentes y en régimen armónico, una función generadora
escalar $\psi$ (una componente del campo) satisface la ecuación de Helmholtz:

$$
\nabla^2 \psi \;+\; \varepsilon_{\text{reg}}\, k^2\, \psi \;=\; 0
\tag{I.1}
$$

donde $\nabla^2 = \partial_{xx}+\partial_{yy}+\partial_{zz}$ es el laplaciano, y
$\varepsilon_{\text{reg}}$ vale $\varepsilon_r$ en la losa y $1$ en el vacío.

## I.2 — Separación de variables

**Se expone.** Suponemos solución en producto de funciones de una variable:

$$
\psi(x,y,z) \;=\; X(x)\,Y(y)\,Z(z)
\tag{I.2}
$$

Sustituyendo (I.2) en (I.1) y dividiendo por $\psi$, cada término depende de una sola variable,
luego cada uno es constante. Llamando $-k_x^2,-k_y^2,-k_z^2$ a esas constantes
($X''/X=-k_x^2$, etc.), la ecuación se vuelve algebraica —la **relación de dispersión**:

$$
k_x^2 + k_y^2 + k_z^2 \;=\; \varepsilon_{\text{reg}}\, k^2
\tag{I.3}
$$

## I.3 — Problema transversal: las paredes laterales fijan $k_m$

**Se expone qué se impone.** Las paredes $x=0,a$ y $y=0,b$ son PEC: el campo tangencial se
anula. Eso obliga a que $X,Y$ sean senos que se anulan en las paredes,

$$
X(x)=\sin\!\frac{m\pi x}{a}\;\Rightarrow\; k_x=\frac{m\pi}{a},
\qquad
Y(y)=\sin\!\frac{n\pi y}{b}\;\Rightarrow\; k_y=\frac{n\pi}{b}.
$$

Definimos el **número de onda transversal** $k_m$ (el "corte" de la sección $a\times b$):

$$
k_m^2 \;=\; k_x^2 + k_y^2 \;=\; \left(\frac{m\pi}{a}\right)^{\!2} + \left(\frac{n\pi}{b}\right)^{\!2}
\tag{I.4}
$$

**Punto clave:** $k_m$ es el **mismo** en las dos regiones, porque en la interfaz $z=h$ los
campos deben casar para todo $(x,y)$; solo cambia la parte en $z$.

## I.4 — Problema longitudinal: $k_z$ en cada región

**Se expone.** Despejando $k_z$ de (I.3) con (I.4) en cada región (cambia
$\varepsilon_{\text{reg}}$, no cambia $k_m$):

$$
k_{z1}^2 = \varepsilon_r\,k^2 - k_m^2 \quad (\text{dieléctrico},\,0\le z\le h),
\qquad
k_{z2}^2 = k^2 - k_m^2 \quad (\text{vacío},\,h\le z\le c).
\tag{I.5}
$$

En cada región $Z''+k_z^2 Z=0$. Las paredes PEC de arriba y abajo fijan la forma de $Z$:

$$
Z_1(0)=0 \;\Rightarrow\; Z_1(z)=A\,\sin(k_{z1}\,z),
\qquad
Z_2(c)=0 \;\Rightarrow\; Z_2(z)=B\,\sin\!\big(k_{z2}(c-z)\big),
$$

con $A,B$ las amplitudes (constantes desconocidas de cada región).

## I.5 — Condiciones de acoplamiento en la interfaz $z=h$

**Se exponen.** Con $\mu$ igual en ambos medios, la continuidad de los campos tangenciales
equivale a que $\psi$ y su derivada en $z$ sean continuas:

$$
Z_1(h)=Z_2(h),
\tag{I.6a}
$$
$$
Z_1'(h)=Z_2'(h).
\tag{I.6b}
$$

Evaluando cada lado (la prima es $d/dz$; ojo al signo al derivar $\sin(k_{z2}(c-z))$):

$$
\begin{aligned}
Z_1(h)&=A\sin(k_{z1}h), & Z_1'(h)&=A\,k_{z1}\cos(k_{z1}h),\\
Z_2(h)&=B\sin\!\big(k_{z2}(c-h)\big), & Z_2'(h)&=-B\,k_{z2}\cos\!\big(k_{z2}(c-h)\big).
\end{aligned}
$$

Con lo que (I.6a) y (I.6b) quedan:

$$
A\sin(k_{z1}h)=B\sin\!\big(k_{z2}(c-h)\big),
\qquad
A\,k_{z1}\cos(k_{z1}h)=-B\,k_{z2}\cos\!\big(k_{z2}(c-h)\big).
$$

## I.6 — Ecuación de resonancia (el resultado)

**Se obtiene eliminando las amplitudes.** Dividiendo la segunda ecuación entre la primera,
$A$ y $B$ desaparecen y, usando $\cot=\cos/\sin$, se llega a la ecuación de resonancia
trascendente:

$$
\boxed{\;
k_{z1}\,\cot\!\big(k_{z1}\,h\big) \;+\; k_{z2}\,\cot\!\big(k_{z2}\,(c-h)\big) \;=\; 0
\;}
\tag{I.7}
$$

$$
k_{z1}=\sqrt{\varepsilon_r\,k^2-k_m^2},\qquad
k_{z2}=\sqrt{k^2-k_m^2},\qquad
k_m^2=\left(\tfrac{m\pi}{a}\right)^2+\left(\tfrac{n\pi}{b}\right)^2.
$$

La **única** incógnita es $k=k_c$; el resto son datos. Es exactamente la ecuación de
`getKcFromHalfFilledCavity.m` y de `_resonance_eq` en `HalfFilledCavity.jl`.

## I.7 — Cómo se resuelve (I.7)

- No es despejable ($k$ dentro de raíces y cotangentes) $\Rightarrow$ se resuelve
  **numéricamente**.
- Fijado $(m,n)$, la ecuación tiene **infinitas raíces** en $k_c$ (una por modo longitudinal),
  porque $\cot$ es periódica.
- A cada raíz se le da una **semilla** cercana. En `mainHalfFilled.m`:

  | modo $(m,n)$ | semillas de $k_c$ |
  |---|---|
  | $(1,0)$ | 300, 600, 800 |
  | $(2,0)$ | 540, 760 |
  | $(3,0)$ | 750 |

- MATLAB usa `vpasolve`; Julia usa bisección (`analytic_kc`).
- Si $k^2<k_m^2/\varepsilon$, la raíz $k_z$ es imaginaria (modo evanescente): $k_z\cot(k_zL)$
  sigue siendo real ($=\kappa\coth(\kappa L)$), por eso el código trabaja en complejos y toma
  la parte real.

Estas raíces $k_c$ son la **solución exacta de referencia**.

---

# PARTE II — Solución por FEM

Idea: cuando la geometría no es una caja simple, no hay fórmula. El FEM resuelve Maxwell de
forma aproximada partiendo el dominio en elementos y convirtiéndolo en álgebra lineal. Aquí se
aplica a la **misma** cavidad para poder comparar con la Parte I.

## II.1 — Forma fuerte (la EDP original)

**Se expone.** Combinando las ecuaciones de Maxwell se obtiene la ecuación de onda vectorial
(doble rotacional) para el campo eléctrico:

$$
\nabla\times\!\left(\mu_r^{-1}\,\nabla\times\mathbf{E}\right) \;=\; k_c^2\,\varepsilon_r\,\mathbf{E}
\tag{II.1}
$$

donde $\mu_r$ es la permeabilidad relativa ($=1$ aquí), $\varepsilon_r$ la permitividad
($2$ en la losa, $1$ en el vacío) y $k_c^2$ el autovalor buscado. Condición PEC en las paredes
(la componente tangencial se anula, con $\hat{\mathbf n}$ la normal a la pared):

$$
\hat{\mathbf n}\times\mathbf{E}=0 \quad\text{en la frontera.}
\tag{II.2}
$$

## II.2 — Forma débil (variacional)

**Se expone por qué.** El FEM no resuelve (II.1) directamente, sino su forma débil:
multiplicamos por una función de test $\mathbf v$ e integramos, bajando el orden de las
derivadas con la identidad de Green:

$$
\underbrace{\int_\Omega \mu_r^{-1}\,(\nabla\times\mathbf E)\cdot(\nabla\times\mathbf v)\,dV}_{a(\mathbf E,\mathbf v)\;=\;\text{RIGIDEZ}}
\;=\;
k_c^2\underbrace{\int_\Omega \varepsilon_r\,\mathbf E\cdot\mathbf v\,dV}_{b(\mathbf E,\mathbf v)\;=\;\text{MASA}}
\tag{II.3}
$$

$\mathbf v$ es una función de test de la misma clase que $\mathbf E$ (tangencial nula en el
borde), $\Omega$ el dominio y $dV$ el diferencial de volumen. La forma $a$ mide "energía de
rotacional"; la forma $b$, "energía del campo" pesada por $\varepsilon_r$.

## II.3 — Discretización con las bases de Nédélec

**Se expone la aproximación.** El campo se aproxima como combinación de funciones de base
conocidas ($54$ por hexaedro, orden 2):

$$
\mathbf E \;\approx\; \sum_{j=1}^{54} x_j\,\mathbf N_j
\tag{II.4}
$$

$\mathbf N_j$ es la $j$-ésima base vectorial de Nédélec (curl-conforme) y $x_j$ su coeficiente
(incógnita). Las $\mathbf N_j$ se generan resolviendo una vez $G\,W=I$ (`get_coefficients`).

## II.4 — Problema de autovalores generalizado

**Se expone el resultado.** Sustituyendo (II.4) en (II.3) y tomando $\mathbf v=\mathbf N_i$, las
integrales se vuelven matrices:

$$
\boxed{\;K\,\mathbf x \;=\; k_c^2\,M\,\mathbf x\;}
\tag{II.5}
$$

$\mathbf x$ es el vector de coeficientes (el modo discreto), $K_{ij}=a(\mathbf N_i,\mathbf N_j)$
la matriz global de rigidez, $M_{ij}=b(\mathbf N_i,\mathbf N_j)$ la de masa, y $k_c^2$ el
autovalor. Es **generalizado** (aparece $M$, no la identidad) porque las bases no son
ortonormales: $\int \mathbf N_i\cdot\mathbf N_j\neq\delta_{ij}$.

## II.5 — Matrices elementales (de dónde salen $K$ y $M$)

**Se exponen.** $K$ y $M$ se montan sumando contribuciones de cada elemento. En un elemento se
calculan dos matrices $54\times54$:

$$
M^e_{ij}=\int_{\Omega_e}\varepsilon_r\,\mathbf N_i\cdot\mathbf N_j\,dV,
\qquad
K^e_{ij}=\int_{\Omega_e}\mu_r^{-1}\,(\nabla\times\mathbf N_i)\cdot(\nabla\times\mathbf N_j)\,dV.
\tag{II.6}
$$

Las integrales se hacen por **cuadratura de Gauss** (no analíticamente):

$$
\int_{\Omega_e} f\,dV \;\approx\; \sum_{g} w_g\,f(\boldsymbol\xi_g)\,\lvert\det J\rvert,
$$

con $w_g$ el peso de Gauss del punto $g$ (27 puntos, orden 4), $\boldsymbol\xi_g$ el punto en el
elemento de referencia $[0,1]^3$, $J$ el Jacobiano del mapeo referencia $\to$ físico y
$\lvert\det J\rvert$ el factor de cambio de volumen. En cada punto se aplica el mapeo de Piola
(covariante para $\mathbf N$, contravariante para $\nabla\times\mathbf N$). Código:
`mass_matrix_element`, `stiff_matrix_element`. **$\varepsilon_r$ entra solo en $M$.**

## II.6 — Ensamblaje global

**Se expone la operación.** Cada DOF local tiene un número global; la entrada local $(i,j)$ de
un elemento se **suma** en la posición global:

$$
K\big[g_i,\,g_j\big] \mathrel{+}= K^e_{ij},
\qquad
M\big[g_i,\,g_j\big] \mathrel{+}= M^e_{ij},
\tag{II.7}
$$

donde $g_i=\texttt{dof\_map.connectivity}[e,i]$ es el número global del DOF local $i$. Como
aristas y caras se comparten, varios elementos suman a la misma casilla $\Rightarrow$ así se
"cosen" los elementos y se impone la continuidad tangencial. Se usa formato COO y `sparse()`
suma los duplicados. Código: `assemble_global_matrices` / `getGlobalMatrices.m`.

## II.7 — Condiciones PEC y resolución

**Se expone.** Los DOFs sobre las paredes (PEC) se eliminan; queda el sistema reducido a los
DOFs libres, $K_f=K[\text{free},\text{free}]$, $M_f=M[\text{free},\text{free}]$. El autovalor
es el **cociente de Rayleigh**:

$$
k_c^2 \;=\; \frac{\mathbf x^{\mathsf T} K\,\mathbf x}{\mathbf x^{\mathsf T} M\,\mathbf x}.
\tag{II.8}
$$

Cómo se resuelve numéricamente:

- $n_\text{free}\le 2000$: `eigvals` denso (todo el espectro) + filtro de espurios.
- $n_\text{free}>2000$: `eigs` ARPACK shift-invert con $\sigma$ centrado en la banda física
  (`explicittransform=:none`).

Código: `_lowest_physical_eigenvalues` (`HalfFilledCavity.jl`) / `eigs` (`mainHalfFilled.m`).

## II.8 — Modos espurios (autovalores nulos)

**Se expone el fenómeno.** El operador $\nabla\times(\nabla\times\cdot)$ tiene un núcleo enorme:
todo campo gradiente $\mathbf E=\nabla\phi$ cumple $\nabla\times\mathbf E=0$, luego
$K\mathbf x=0$. Esos modos dan $k_c^2\approx 0$:

$$
\text{espurios: } k_c^2\approx 0 \quad(\text{campos gradiente}),
\qquad
\text{físicos: } k_c^2\sim 10^{5}.
$$

Se filtran con umbral fijo ($k_c^2>1$) y se toman los más bajos físicos.

## II.9 — Validación: convergencia

**Se expone qué se espera.** Para un elemento de orden $p=2$, el error del autovalor decae

$$
\big|\,\text{error}(k_c^2)\,\big| \;=\; O\!\big(h^{2p}\big) \;=\; O\!\big(h^{4}\big),
\tag{II.9}
$$

con $h$ el tamaño característico del elemento. (El autovalor converge al doble de velocidad que
el campo, $O(h^2)$, por superconvergencia del cociente de Rayleigh.) La curva log–log de
pendiente $\approx 4$ valida la implementación (paper Fig. 7, `test_10`).

---

# El puente entre las dos soluciones

$$
\text{error relativo} \;=\; \frac{\big|\,k_{c,\text{FEM}}^2 - k_{c,\text{analítico}}^2\,\big|}{\big|\,k_{c,\text{analítico}}^2\,\big|}
\;\xrightarrow[\ h\to 0\ ]{}\; 0 \quad\text{como } O(h^4).
$$

- **Parte I (analítica):** raíces $k_c$ de (I.7) $\Rightarrow$ verdad exacta.
- **Parte II (FEM):** autovalores $k_c^2$ de (II.5) $\Rightarrow$ aproximación.

**En una frase:** la solución analítica resuelve Maxwell por separación de variables e impone
las paredes y la continuidad en la interfaz para llegar a la ecuación trascendente
$k_{z1}\cot(k_{z1}h)+k_{z2}\cot(k_{z2}(c-h))=0$, cuyas raíces $k_c$ son exactas; la solución FEM
convierte la misma física en la forma débil, la discretiza con las 54 bases de Nédélec en el
problema $K\mathbf x=k_c^2 M\mathbf x$, lo resuelve filtrando los modos espurios, y se valida
comprobando que sus autovalores convergen $O(h^4)$ a los $k_c$ analíticos.
