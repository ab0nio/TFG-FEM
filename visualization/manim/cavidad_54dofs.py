import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from manim import *

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# ─── Geometría real del benchmark (HalfFilledCavity.jl) ──────────────────────
# a=10mm, b=1mm, c=10mm, losa dieléctrica εr=2 en z ∈ [0, h=5mm].
# Malla del test 4: nx=4, ny=1, nz=4. Escala de escena: 1 mm -> 0.48 unidades.
A_MM, B_MM, C_MM, H_MM = 10.0, 1.0, 10.0, 5.0
S = 0.48
NX, NY, NZ = 4, 1, 4

# Semilados en unidades de escena (caja centrada en el origen, z hacia arriba)
HX, HY, HZ = A_MM * S / 2, B_MM * S / 2, C_MM * S / 2

# Topología del hexaedro (convención ReferenceHex.jl, índices base 0)
EDGES = [
    (0, 1), (1, 2), (3, 2), (0, 3),
    (4, 5), (5, 6), (7, 6), (4, 7),
    (0, 4), (1, 5), (2, 6), (3, 7),
]
UNIT_NODES = np.array([
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
], dtype=float)

REF_SIDE = 2.6  # lado del hexaedro de referencia en pantalla


def ref(p):
    """[0,1]^3 -> escena, centrado en el origen."""
    return (np.asarray(p, dtype=float) - 0.5) * REF_SIDE


def wireframe(corners, color, stroke_width=2.5, opacity=1.0):
    return VGroup(*[
        Line(corners[i], corners[j], color=color,
             stroke_width=stroke_width, stroke_opacity=opacity)
        for i, j in EDGES
    ])


def box_corners(x0, x1, y0, y1, z0, z1):
    return [np.array([x0 + (x1 - x0) * n[0],
                      y0 + (y1 - y0) * n[1],
                      z0 + (z1 - z0) * n[2]]) for n in UNIT_NODES]


def glyph(center, direction, color, length=0.55):
    """Flecha barata para Cairo: segmento + punto en la punta.
    (Arrow3D es demasiado costoso en masa, ver continuidad_tangencial.py)"""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d) * length
    p0, p1 = center - d / 2, center + d / 2
    return VGroup(
        Line(p0, p1, color=color, stroke_width=4),
        Dot(p1, radius=0.05, color=color),
    )


def heatmap_image(mat, rgb, height=4.0):
    """Matriz 54x54 como una única textura (mucho más barato que 2916 Squares)."""
    v = np.abs(mat) / np.abs(mat).max()
    img = np.zeros((*mat.shape, 3), dtype=np.uint8)
    for k in range(3):
        img[..., k] = (v * rgb[k]).astype(np.uint8)
    im = ImageMobject(img)
    im.set_resampling_algorithm(RESAMPLING_ALGORITHMS["nearest"])
    im.height = height
    return im


class CavidadSemillena54DOFs(ThreeDScene):
    """
    Pipeline completo del TFG en una escena:
      1. Cavidad semillena real (a=10mm, b=1mm, c=10mm, εr=2 en z<5mm)
      2. Mallado estructurado 4x1x4 (la malla del test del benchmark)
      3. Zoom a un elemento -> hexaedro de referencia [0,1]^3
      4. Los 54 DOFs en sus posiciones (24 arista + 24 cara + 6 volumen)
      5. Las integrales funcionales de cada tipo (ecs. 6-8 del paper)
      6. Construcción de G (G_ij = g_i(u_j)) y W = G^-1 con las matrices
         REALES exportadas de CoefficientGenerator.jl (data/*.csv)
    """

    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=-70 * DEGREES)

        # ══ 1. Cavidad semillena ══════════════════════════════════════════
        title = Text("Cavidad semillena", font_size=34).to_corner(UL)
        subtitle = Text("a = 10 mm,  b = 1 mm,  c = 10 mm", font_size=22,
                        color=GRAY_B).next_to(title, DOWN, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(title, subtitle)

        cavity = wireframe(box_corners(-HX, HX, -HY, HY, -HZ, HZ), WHITE, 3)
        dielectric = Prism(dimensions=[2 * HX, 2 * HY, HZ],
                           fill_color=ORANGE, fill_opacity=0.3, stroke_width=0)
        dielectric.move_to([0, 0, -HZ / 2])

        legend_die = Text("dieléctrico  εr = 2   (z < h = 5 mm)",
                          font_size=22, color=ORANGE).to_corner(DL)
        legend_vac = Text("vacío  εr = 1", font_size=22,
                          color=BLUE_B).next_to(legend_die, UP, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(legend_die, legend_vac)

        self.play(Create(cavity), Write(title), FadeIn(subtitle), run_time=2)
        self.play(FadeIn(dielectric), FadeIn(legend_die), FadeIn(legend_vac),
                  run_time=1.5)
        self.wait(0.5)

        # ══ 2. Mallado 4x1x4 ══════════════════════════════════════════════
        mesh_label = Text("Malla estructurada 4×1×4  (16 hexaedros)",
                          font_size=24, color=TEAL).to_corner(UR)
        self.add_fixed_in_frame_mobjects(mesh_label)

        xs = np.linspace(-HX, HX, NX + 1)
        ys = np.linspace(-HY, HY, NY + 1)
        zs = np.linspace(-HZ, HZ, NZ + 1)
        mesh = VGroup()
        for y in ys:
            for z in zs:
                mesh.add(Line([-HX, y, z], [HX, y, z]))
        for x in xs:
            for z in zs:
                mesh.add(Line([x, -HY, z], [x, HY, z]))
        for x in xs:
            for y in ys:
                mesh.add(Line([x, y, -HZ], [x, y, HZ]))
        mesh.set_stroke(TEAL, width=1.5, opacity=0.9)

        self.play(FadeIn(mesh_label))
        self.play(Create(mesh, lag_ratio=0.02), run_time=3)
        self.wait(0.5)

        # ══ 3. Zoom a un elemento ═════════════════════════════════════════
        # Elemento con x ∈ [0, a/4], z ∈ [-c/4, 0] (capa superior del dieléctrico)
        elem = wireframe(box_corners(0, HX / 2, -HY, HY, -HZ / 2, 0),
                         YELLOW, 5)
        self.play(Create(elem), run_time=1)
        self.play(
            cavity.animate.set_stroke(opacity=0.15),
            mesh.animate.set_stroke(opacity=0.1),
            FadeOut(dielectric),
        )
        self.wait(0.3)

        ref_cube = wireframe([ref(n) for n in UNIT_NODES], WHITE, 3)
        ref_title = Text("Hexaedro de referencia  [0,1]³", font_size=30).to_corner(UL)
        self.add_fixed_in_frame_mobjects(ref_title)
        ref_title.set_opacity(0)

        self.play(
            ReplacementTransform(elem, ref_cube),
            FadeOut(cavity), FadeOut(mesh),
            FadeOut(title), FadeOut(subtitle),
            FadeOut(legend_die), FadeOut(legend_vac), FadeOut(mesh_label),
            ref_title.animate.set_opacity(1),
            run_time=2,
        )
        self.wait(0.5)

        # ══ 4. Los 54 DOFs en sus posiciones ══════════════════════════════
        # Aristas: 2 DOFs por arista (q = 1-t y q = t), glifos en t=0.3 y t=0.7
        edge_glyphs = VGroup()
        for a, b in EDGES:
            p1, p2 = ref(UNIT_NODES[a]), ref(UNIT_NODES[b])
            tau = p2 - p1
            for t in (0.3, 0.7):
                edge_glyphs.add(glyph(p1 + tau * t, tau, RED, length=0.5))

        # Caras: (eje fijo, valor, eje u, eje v) en el orden de DOFs.jl
        face_defs = [
            (2, 0.0, 0, 1), (2, 1.0, 0, 1),   # z=0, z=1: varían x,y
            (1, 0.0, 0, 2),                    # y=0: varían x,z
            (0, 1.0, 1, 2),                    # x=1: varían y,z
            (1, 1.0, 0, 2),                    # y=1: varían x,z
            (0, 0.0, 1, 2),                    # x=0: varían y,z
        ]
        face_glyphs = VGroup()
        for fixed_ax, fixed_val, u_ax, v_ax in face_defs:
            # 4 DOFs por cara: 2 con q a lo largo del eje u, 2 a lo largo del eje v
            for ax, positions in ((u_ax, [(1/3, 1/2), (2/3, 1/2)]),
                                  (v_ax, [(1/2, 1/3), (1/2, 2/3)])):
                for pu, pv in positions:
                    p = np.full(3, 0.5)
                    p[fixed_ax], p[u_ax], p[v_ax] = fixed_val, pu, pv
                    d = np.zeros(3)
                    d[ax] = 1.0
                    face_glyphs.add(glyph(ref(p), d, BLUE, length=0.5))

        # Volumen: 6 DOFs, q = (1-ξ)x̂, ξx̂, (1-η)ŷ, ηŷ, (1-ζ)ẑ, ζẑ
        vol_glyphs = VGroup()
        for ax in range(3):
            for t in (1/3, 2/3):
                p = np.full(3, 0.5)
                p[ax] = t
                d = np.zeros(3)
                d[ax] = 1.0
                vol_glyphs.add(glyph(ref(p), d, GREEN, length=0.5))

        lbl_edge = Text("Aristas: 12 × 2 = 24", font_size=24, color=RED).to_corner(UR)
        lbl_face = Text("Caras: 6 × 4 = 24", font_size=24, color=BLUE)
        lbl_vol = Text("Volumen: 1 × 6 = 6", font_size=24, color=GREEN)
        lbl_face.next_to(lbl_edge, DOWN, aligned_edge=RIGHT)
        lbl_vol.next_to(lbl_face, DOWN, aligned_edge=RIGHT)
        self.add_fixed_in_frame_mobjects(lbl_edge, lbl_face, lbl_vol)
        VGroup(lbl_edge, lbl_face, lbl_vol).set_opacity(0)

        self.play(lbl_edge.animate.set_opacity(1),
                  Create(edge_glyphs, lag_ratio=0.05), run_time=2.5)
        self.play(lbl_face.animate.set_opacity(1),
                  Create(face_glyphs, lag_ratio=0.05), run_time=2.5)
        self.play(lbl_vol.animate.set_opacity(1),
                  Create(vol_glyphs, lag_ratio=0.1), run_time=1.5)

        total = MathTex(r"24 + 24 + 6 = 54\ \text{DOFs}", font_size=38).to_edge(DOWN)
        self.add_fixed_in_frame_mobjects(total)
        total.set_opacity(0)
        self.play(total.animate.set_opacity(1))
        self.wait(1)
        self.play(total.animate.set_opacity(0))

        # ══ 5. Integrales funcionales (ecs. 6-8 del paper) ════════════════
        eq_edge = MathTex(
            r"g_i(\mathbf{u}) = \int_{e}(\mathbf{u}\cdot\hat{\tau})\,q\,dt,"
            r"\qquad q \in \{1-t,\; t\}",
            font_size=34, color=RED,
        ).to_edge(DOWN)
        eq_face = MathTex(
            r"g_i(\mathbf{u}) = \iint_{f}(\hat{n}\times\mathbf{u})\cdot\mathbf{q}\,dS,"
            r"\qquad \text{4 } \mathbf{q} \text{ tangentes por cara}",
            font_size=34, color=BLUE,
        ).to_edge(DOWN)
        eq_vol = MathTex(
            r"g_i(\mathbf{u}) = \iiint_{V}\mathbf{u}\cdot\mathbf{q}\,dV,"
            r"\qquad \mathbf{q} \in \{(1-\xi)\hat{x},\,\xi\hat{x},\,\ldots\}",
            font_size=34, color=GREEN,
        ).to_edge(DOWN)

        for eq, grp, others in (
            (eq_edge, edge_glyphs, VGroup(face_glyphs, vol_glyphs)),
            (eq_face, face_glyphs, VGroup(edge_glyphs, vol_glyphs)),
            (eq_vol, vol_glyphs, VGroup(edge_glyphs, face_glyphs)),
        ):
            self.add_fixed_in_frame_mobjects(eq)
            eq.set_opacity(0)
            self.play(eq.animate.set_opacity(1),
                      grp.animate.set_stroke(opacity=1).set_fill(opacity=1),
                      others.animate.set_stroke(opacity=0.12).set_fill(opacity=0.12))
            self.wait(1.5)
            self.play(eq.animate.set_opacity(0))
        self.play(VGroup(edge_glyphs, face_glyphs, vol_glyphs)
                  .animate.set_stroke(opacity=1).set_fill(opacity=1))
        self.wait(0.3)

        # ══ 6. De los DOFs a la matriz de coeficientes ════════════════════
        self.play(
            FadeOut(ref_cube), FadeOut(edge_glyphs), FadeOut(face_glyphs),
            FadeOut(vol_glyphs), FadeOut(ref_title),
            FadeOut(lbl_edge), FadeOut(lbl_face), FadeOut(lbl_vol),
        )

        G = np.loadtxt(os.path.join(DATA_DIR, "G_matrix.csv"), delimiter=",")
        W = np.loadtxt(os.path.join(DATA_DIR, "W_matrix.csv"), delimiter=",")

        step_title = Text("Obtención de las funciones base (Fig. 2 del paper)",
                          font_size=28).to_edge(UP)
        eq_G = MathTex(r"G_{ij} = g_i(u_j)", font_size=40,
                       substrings_to_isolate=["g_i", "u_j"])
        eq_G.set_color_by_tex("g_i", RED)
        eq_G.set_color_by_tex("u_j", YELLOW)
        note_G = Text("54 funcionales × 54 monomios de H₂",
                      font_size=22, color=GRAY_B)

        g_img = heatmap_image(G, (80, 150, 255))
        g_cap = Text("G", font_size=28, color=BLUE_B)

        self.add_fixed_in_frame_mobjects(step_title, eq_G, note_G, g_img, g_cap)
        eq_G.next_to(step_title, DOWN, buff=0.4)
        note_G.next_to(eq_G, DOWN, buff=0.2)
        g_img.next_to(note_G, DOWN, buff=0.35)
        g_img.shift(LEFT * 2.8)
        g_cap.next_to(g_img, DOWN, buff=0.15)

        for m in (step_title, eq_G, note_G, g_cap):
            m.set_opacity(0)
        g_img.set_opacity(0)

        self.play(step_title.animate.set_opacity(1), eq_G.animate.set_opacity(1),
                  note_G.animate.set_opacity(1))
        self.play(g_img.animate.set_opacity(1), g_cap.animate.set_opacity(1),
                  run_time=1.5)
        self.wait(0.8)

        eq_W = MathTex(r"G\,W = I_{54} \;\Rightarrow\; W = G^{-1}", font_size=40)
        w_img = heatmap_image(W, (90, 230, 140))
        w_cap = Text("W  (coeficientes de las 54 bases)", font_size=24, color=GREEN_B)

        self.add_fixed_in_frame_mobjects(eq_W, w_img, w_cap)
        w_img.move_to(g_img.get_center() + RIGHT * 5.6)
        w_cap.next_to(w_img, DOWN, buff=0.15)
        arrow = Arrow(g_img.get_right(), w_img.get_left(), buff=0.25, color=WHITE)
        self.add_fixed_in_frame_mobjects(arrow)
        eq_W.move_to(eq_G)  # sustituye a G_ij en la cabecera, sin tapar los heatmaps

        for m in (eq_W, w_cap, arrow):
            m.set_opacity(0)
        w_img.set_opacity(0)
        self.play(eq_G.animate.set_opacity(0), note_G.animate.set_opacity(0),
                  eq_W.animate.set_opacity(1), arrow.animate.set_opacity(1))
        self.play(w_img.animate.set_opacity(1), w_cap.animate.set_opacity(1),
                  run_time=1.5)
        self.wait(1)

        # Bloque 5x5 real de W: los coeficientes son ENTEROS (Tabla II)
        block = [[str(int(round(W[i, j]))) for j in range(5)] for i in range(5)]
        mat_tex = Matrix(block, h_buff=1.1).scale(0.75)
        block_cap = Text("Bloque real de W: coeficientes enteros (Tabla II del paper)",
                         font_size=24, color=YELLOW)
        final_grp = VGroup(mat_tex, block_cap).arrange(DOWN, buff=0.45)
        self.add_fixed_in_frame_mobjects(final_grp)
        final_grp.set_opacity(0)

        self.play(
            g_img.animate.set_opacity(0), w_img.animate.set_opacity(0),
            *[m.animate.set_opacity(0) for m in
              (eq_G, eq_W, g_cap, w_cap, arrow, step_title)],
        )
        self.play(final_grp.animate.set_opacity(1), run_time=1.5)
        self.wait(2.5)
