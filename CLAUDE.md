# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TFG (Trabajo de Fin de Grado) implementing the Finite Element Method, specifically second-order Nédélec curl-conforming hexahedral elements. The repository has two independent implementations:

- **`NedelecHex/`** — Julia package (primary, active codebase)
- **`matlab/`** — MATLAB implementation (reference/legacy, used for numerical verification against a half-filled cavity benchmark)

There is also a **`visualization/`** directory with Python/Manim scripts for animations.

## Julia Package: NedelecHex

### Commands

```bash
# Enter Julia REPL in the package environment
cd NedelecHex
julia --project=.

# Run all tests
julia --project=. test/runtests.jl

# Run all tests in verbose mode: each block prints what it checks, why the
# expected value is correct, and the measured quantities
NEDELEC_TEST_VERBOSE=1 julia --project=. test/runtests.jl

# Cross-comparison with MATLAB output (requires running mainDemo.m first)
julia --project=. demo.jl

# Evaluate basis functions at the reference element center
julia --project=. demo_bases.jl

# Convergence study on MATLAB-exported meshes (requires CSVs from exportMesh.m)
using NedelecHex
kc = NedelecHex.HalfFilledCavity.analytic_eigenvalues()
errors = NedelecHex.HalfFilledCavity.convergence_study_matlab(kc; mesh_dir="path/to/csvs")
```

Inside the Julia REPL:
```julia
using Pkg; Pkg.test()  # run all tests
```

### Architecture

The module load order in `src/NedelecHex.jl` reflects the dependency graph:

1. **`Geometry/ReferenceHex.jl`** — reference hexahedral element geometry (no deps). Defines `NODES` (8 corners of [0,1]³), `EDGES` (12 edge pairs), `FACES` (6 face 4-tuples), `FACE_NORMALS`. **All other modules assume this node ordering.**
2. **`Geometry/Mesh.jl`** — generates structured hex meshes via Gmsh transfinite. Exports `HexMesh` and `box_mesh`. **Critical:** Gmsh's transfinite hex8 node output uses a rotated ordering; `box_mesh` applies the permutation `[3,7,8,4,2,6,5,1]` after building the connectivity to match `ReferenceHex.NODES`. Without this, the Jacobian is a rotation matrix instead of diagonal, breaking DOF assembly.
3. **`Math/GaussQuadrature.jl`** — Gauss–Legendre quadrature points and weights on [0,1]. `gauss_hex(order)` returns a tensor-product 3D rule.
4. **`Math/CoordinateTransforms.jl`** — isoparametric mapping between reference [0,1]³ and physical coordinates.
5. **`Bases/Monomials.jl`** — the 54 vectorial monomials spanning H₂ and their analytical curl vectors. `build_h2_monomials()` and `build_h2_curl_monomials()`.
6. **`Bases/DOFs.jl`** — the 54 DOF functionals (edge, face, volume) as defined in Table I of the paper. `Q_EDGE`, `Q_FACE`, `Q_VOL` test vectors; `evaluateEdgeDOF`, `evaluateFaceDOF`, `evaluateVolumeDOF`.
7. **`Bases/CoefficientGenerator.jl`** — builds the 54×54 DOF matrix G and solves W = G\I. `get_coefficients()` returns `(W, G)`. `run_all_tests()` verifies W against Table II of the paper.
8. **`Bases/BasisEval.jl`** — evaluates basis functions and their curls at a reference point. `eval_bases_ref`, `eval_curl_bases_ref`, `compute_jacobian`, `piola_transform`, `curl_piola_transform`.
9. **`Assembly/GlobalDOFs.jl`** — builds the global DOF connectivity table (n_elem × 54). Handles edge DOF orientation (swap when local edge runs opposite to canonical min→max direction). Face DOFs identified by sorted 4-node key. Volume DOFs are element-local. `build_dof_map`, `dirichlet_dofs`.
10. **`Assembly/Assembly.jl`** — assembles sparse global stiffness K and mass M via COO format. Integration order 4 (3-point Gauss rule, 27 points/element). `assemble_global_matrices`, `free_dofs`.
11. **`Benchmark/HalfFilledCavity.jl`** — half-filled cavity benchmark: analytic eigenvalues via bisection on the transcendental resonance equation (`analytic_kc`, `analytic_eigenvalues`), per-element permittivity (`_element_epsilon`), and `convergence_study_matlab` (runs the FEM benchmark on meshes exported from MATLAB with `exportMesh.m`, for a 1:1 comparison with `mainHalfFilled.m`). The self-contained Julia equivalent of the `mainHalfFilled.m` pipeline (mesh → assembly → eigenvalues vs analytic) lives in `test/test_10_half_filled_cavity.jl`.

### Test Suite

`test/runtests.jl` includes 10 files ordered by dependency (each level relies only on the previous ones). With `NEDELEC_TEST_VERBOSE=1`, every block prints what it checks and why the expected value is correct, via the `explica()` helper defined in `runtests.jl`.

| File | What it checks |
|---|---|
| `test_01_quadrature.jl` | Gauss rules 1D/2D/3D: weights sum to the measure, points in range, exact integration of ξ²η²ζ² |
| `test_02_monomials.jl` | 54 monomials: single nonzero component per block, Nédélec degree constraint via n-th finite differences, analytic curls vs centered-difference curls |
| `test_03_dofs.jl` | Edge/face/volume DOF functionals vs hand-computed integrals of nontrivial polynomial fields |
| `test_04_coefficients.jl` | G unisolvent (rank 54, moderate cond), structural zeros by axis compatibility, G·W = I, W integer ≤ 216 and equal to the full 54×54 Table II (`test/data/tabla2_paper.txt`) |
| `test_05_reference_element.jl` | Jacobian on known geometries, dual property gᵢ(Nⱼ)=δᵢⱼ end-to-end, null modes of Kᵉ have curl ≡ 0, Kᵉ/Mᵉ symmetry/definiteness, Table III condition numbers (19 and 35) |
| `test_06_gmsh_node_order.jl` | The `[3,7,8,4,2,6,5,1]` permutation verified against raw Gmsh output; diagonal Jacobian on every element of a 3×2×4 mesh |
| `test_07_assembly.jl` | Single-element and 2×2×2 assembly: combinatorial DOF counts, symmetry, SPD mass, K kernel dimension 26 = 27−1 (gradients of the triquadratic nodal basis) |
| `test_08_two_elements.jl` | Shared-face DOF identification (96 = 54+54−12), COO accumulation verified entry by entry, consistent edge orientation across the shared face |
| `test_09_mms.jl` | MMS: exact recovery of a manufactured solution inside H₂ (error < 1e-9) and log-log convergence slope ≈ 2 for an oblique plane wave |
| `test_10_half_filled_cavity.jl` | Analytic resonance roots (residual + sign change), FEM eigenvalues on a coarse mesh with spurious-mode filtering, O(h⁴) eigenvalue convergence (paper Fig. 7) |

### Key Implementation Details

**Piola transforms** (`BasisEval.jl`):

The Jacobian `J` stored internally satisfies `J[i,j] = Σ_k (∂N_k/∂ξ_i)(pt) × x_k[j]`, which equals the transpose of the standard isoparametric Jacobian.

- **H(curl) covariant Piola** (for basis functions): `N_phys = N_ref × J⁻¹`, i.e. `(J \ N_ref')'`. This preserves tangential continuity across elements.
- **Curl contravariant Piola** (for curl of basis functions): `curlN_phys = curlN_ref × J / det(J)`. This preserves the curl of the field.

**Global DOF assembly** (`GlobalDOFs.jl`):

Edge DOFs use the canonical orientation min→max global node. If a local edge runs in the opposite direction, the two DOFs are swapped so that DOF₁ corresponds to weight `(1−t)` and DOF₂ to weight `t` consistently.

Face DOFs are shared between the two elements touching a face, identified by the sorted tuple of 4 global node indices.

**Eigenvalue solver** (`HalfFilledCavity.jl`):

The H(curl) stiffness matrix K has a large null space (irrotational/gradient spurious modes with curl = 0). This makes ARPACK's shift-and-invert unreliable for this problem — it converges to the high-frequency modes regardless of the shift value. The implementation uses `LinearAlgebra.eigvals` (dense solver) with selection by proximity to sigma, which is reliable for the mesh sizes in the benchmark (n_free ≤ ~1000 for nx ≤ 8).

**Mesh node ordering** (`Mesh.jl`):

Gmsh's transfinite hex8 mesh stores nodes in a rotated order relative to the `ReferenceHex.NODES` convention `[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]`. The permutation `[3,7,8,4,2,6,5,1]` (applied after reading Gmsh output) remaps local node `k` to the correct slot. Without this fix, the Jacobian for every element is a pure rotation matrix `[[0,-Δx,0],[0,0,-Δy],[Δz,0,0]]` instead of `diag(Δx,Δy,Δz)`, which scrambles all DOF assembly.

### Mejoras de rendimiento pendientes (backlog)

- **`GlobalDOFs.jl`, `dirichlet_dofs`, línea de caras**: `shared = reduce(.&, [on_plane[gn[ln], :] for ln in face_lnodes])` crea un array intermedio (comprehension con `[...]`) antes de reducirlo. Cambiar a un generador — `reduce(.&, on_plane[gn[ln], :] for ln in face_lnodes)` (sin corchetes) — da el mismo resultado (verificado) y mide menos de la mitad de bytes asignados (1280 → 544 en una prueba con `@allocated`). Impacto real en tiempo de ejecución: nulo, ya que `dirichlet_dofs` se ejecuta una sola vez por malla, no dentro del bucle de ensamblaje — es una mejora de calidad/idiomática, no de rendimiento medible.

## MATLAB Implementation

Entry points (run from `matlab/` in MATLAB):
- **`mainHalfFilled.m`** — full pipeline: generates basis functions and solves the half-filled cavity eigenvalue problem (numerical verification benchmark)
- **`mainGenerateCoefficients.m`** — generates only basis function coefficients
- **`mainEvaluateCoefficients.m`** — evaluates basis functions from pre-computed coefficients (Table II)

Package structure uses MATLAB's `+package` namespace convention:
- `+assembly/` — global matrix assembly (mass, stiffness, Dirichlet BCs)
- `+bases/` — polynomial basis and `@SystematicBases` class
- `+geometry/` — mesh OOP classes (`@Mesh`, `@Element`, `@Face`, `@Edge`, `@Node`, etc.)
- `+math/` — Gauss quadrature, polynomial spaces, coordinate transforms
- `+settings/` — material properties and system configuration

## Reference Paper

**"Second-Order Nédélec Curl-Conforming Hexahedral Element for Computational Electromagnetics"**
Adrian Amor-Martin and Luis Emilio Garcia-Castillo, *IEEE Trans. Antennas Propagation*, Vol. 71, No. 1, January 2023.
DOI: 10.1109/TAP.2022.3216554

### Mathematical Foundation

The **reference hexahedron** H = [0,1]³ uses coordinates (ξ, η, ζ). The second-order mixed space **H₂** (eq. 4 in the paper) is a tensor product of polynomial spaces ensuring Nédélec's constraint: each component has lower degree in its own axis direction and higher degree in the others. The basis has **54 functions** (N₁…N₅₄), expressed as integer-coefficient linear combinations of monomials — an advantage for symbolic computation.

### Degrees of Freedom (54 total)

DOFs are defined as functionals acting on the vector field **u**:

1. **Edge DOFs (1–24):** 12 edges × 2 DOFs. Integral of tangential component against scalar q ∈ {1−t, t} where t ∈ {ξ, η, ζ} is the edge parameter (see Table I of the paper).
2. **Face DOFs (25–48):** 6 faces × 4 DOFs. Integral of (n̂×**u**) against vector **q** on the face plane.
3. **Volume DOFs (49–54):** Integral of **u** against vector **q** over the whole element.

The specific q and **q** choices for each DOF are in Table I of the paper and implemented in `DOFs.jl` / `evaluateEdgeDOF`, `evaluateFaceDOF`, `evaluateVolumeDOF` in MATLAB.

### Coefficient Generation Algorithm (Fig. 2 in paper)

Implemented in `CoefficientGenerator.jl` (Julia) and `getCoefficients` (MATLAB):
1. Build 54×54 matrix **G** where G(i,j) = g_i(**u**_j) (DOF functional i applied to monomial j)
2. Solve **GW** = **I₅₄** → coefficient matrix **W** (each column defines one basis function as a linear combination of monomials)

The resulting coefficients are integers (≤ 3 digits), as shown in **Table II** of the paper. The Julia test suite (`run_all_tests()` in `CoefficientGenerator.jl`) verifies against this table.

### Basis Function Evaluation in Real Elements (Fig. 3 in paper)

Implemented in `BasisEval.jl` (Julia) / `getBasesEvaluated` (MATLAB):
1. Map real coordinates (x,y,z) → reference coordinates (ξ,η,ζ) via isoparametric mapping
2. Evaluate monomial row vectors P_ξ, P_η, P_ζ at (ξ,η,ζ)
3. Compute reference-element values: N_{ξ,i} = P_ξ·W_i, etc.
4. Apply **covariant Piola mapping**: **N**_i ← J⁻ᵀ **N**_i (preserves tangential continuity / curl-conformity across elements)

### Verification Strategy

- **Half-filled cavity benchmark:** Generalized eigenvalue problem with known analytic eigenvalues from a transcendental resonance equation; convergence follows O(h²) for the second-order element. This is the main benchmark in both `mainHalfFilled.m` and `HalfFilledCavity.jl`.
- **Method of Manufactured Solutions (MMS):** Implemented in `test/test_09_mms.jl` (Julia only; no MATLAB counterpart). Two levels: exact recovery of a manufactured solution inside H₂ with the forcing term derived from the double-curl PDE (eq. 13), error ≈ machine precision on any mesh; and a genuine convergence study with a transverse plane wave (zero forcing), measured log-log slope ≈ 2.
- **Condition numbers** of preconditioned mass **M^p** and stiffness **K^p** matrices for the reference hexahedron are 19 and 35 respectively (Table III of the paper).

### Julia ↔ Paper Pseudocode Mapping

| Paper | Julia function/file |
|---|---|
| Fig. 2 (`GETCOEFFICIENTS`) | `get_coefficients()` in `CoefficientGenerator.jl` |
| Fig. 3 (`NIEVAL`) | `eval_bases_ref`, `piola_transform` in `BasisEval.jl` |
| Fig. 3, line 11 (Jacobian) | `compute_jacobian` in `BasisEval.jl` |
| DOF eqs. (6)–(8) | `evaluateEdgeDOF`, `evaluateFaceDOF`, `evaluateVolumeDOF` in `DOFs.jl` |
| Bilinear forms c₁, c₂ (eqs. 20–21) | `mass_matrix_element`, `stiff_matrix_element` in `Assembly.jl` |
| Global assembly | `assemble_global_matrices` in `Assembly.jl` |
| PEC boundary conditions | `dirichlet_dofs` in `GlobalDOFs.jl` |
| Resonance equation (half-filled cavity) | `_resonance_eq`, `analytic_kc` in `HalfFilledCavity.jl` |

### MATLAB ↔ Paper Pseudocode Mapping

| Paper | MATLAB function |
|---|---|
| Fig. 2 (`GETCOEFFICIENTS`) | `getCoefficients` in `@SystematicBases` |
| Fig. 3 (`NIEVAL`) | `getBasesEvaluated` |
| Fig. 3, line 11 (Jacobian) | `computeJacobianMatrix` |
| DOF eqs. (6)–(8) | `evaluateEdgeDOF`, `evaluateFaceDOF`, `evaluateVolumeDOF` |
| Bilinear forms c₁, c₂ (eqs. 20–21) | `getMassMatrix`, `getStiffMatrix` |
| Human-readable basis display | `showBases`, `showPolynomials` |

## Visualization

`visualization/manim/` contains Python/Manim scripts (local `.venv` with ManimCE 0.20.1; see `IDEAS.md` there for the backlog of planned scenes and render-cost guidance):
- `cavidad_54dofs.py` — main scene `CavidadSemillena54DOFs`: half-filled cavity → 4×1×4 mesh → zoom to one element → 54 DOFs in position → functional integrals (eqs. 6–8) → G/W coefficient matrices with real Julia data
- `base_dual_tabla2.py` — standalone G → W = G⁻¹ animation with the real matrices
- `continuidad_tangencial.py` — tangential continuity of a face basis function shared by two elements (real `BasisEval.jl` data, max|N_A−N_B| ≈ 8e-14)
- `data/` — CSVs exported from Julia (`G_matrix.csv`, `W_matrix.csv`, `face_continuity.csv`) consumed by the scenes

Run with: `cd visualization/manim && .venv/bin/manim -pql <script>.py <SceneName>`

**Performance (Cairo renderer, laptop without dGPU):** prefer flat VMobjects (`Line`, `Dot`, `Polygon`) even inside `ThreeDScene`; render heatmaps as a single `ImageMobject` texture instead of thousands of `Square`s. Avoid `Surface`, mass `Arrow3D`, and long `ambient_camera_rotation` — every frame re-projects the whole scene.
