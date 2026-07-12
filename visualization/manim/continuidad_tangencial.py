import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from manim import *
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Topología del cubo unitario (para dibujar los dos hexaedros)
EDGES = [
    (0, 1), (1, 2), (3, 2), (0, 3),
    (4, 5), (5, 6), (7, 6), (4, 7),
    (0, 4), (1, 5), (2, 6), (3, 7),
]
UNIT_NODES = [
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
]

SCALE = 1.6
ARROW_SCALE = 0.18


def cube_edges(offset, color=GRAY_B):
    pts = [
        (np.array(UNIT_NODES[i]) + np.array(offset) - np.array([1.0, 0.5, 0.5])) * SCALE
        for i in range(8)
    ]
    return VGroup(*[
        Line3D(pts[a], pts[b], color=color, thickness=0.01) for a, b in EDGES
    ])


class ContinuidadTangencial(ThreeDScene):
    """
    Idea 3: evalúa una función base real de cara (DOF de la cara compartida
    x=1 entre dos hexaedros A=[0,1]^3 y B=[1,2]x[0,1]x[0,1]) desde AMBOS
    elementos usando BasisEval.eval_bases_ref + piola_transform, y comprueba
    que los vectores coinciden en la cara compartida (datos reales
    exportados por Julia en data/face_continuity.csv, max|N_A-N_B| ~ 1e-14).
    """

    def construct(self):
        data = np.loadtxt(os.path.join(DATA_DIR, "face_continuity.csv"), delimiter=",")

        self.set_camera_orientation(phi=70 * DEGREES, theta=35 * DEGREES)

        title = Text("Continuidad de la función base en la cara compartida", font_size=28)
        title.to_corner(UL)
        self.add_fixed_in_frame_mobjects(title)

        # ── Dos elementos separados ────────────────────────────────────
        cube_A = cube_edges((0, 0, 0), color=BLUE_D)
        cube_B = cube_edges((2.2, 0, 0), color=ORANGE)  # separados en x

        label_A = Text("Elemento A", font_size=24, color=BLUE_D)
        label_B = Text("Elemento B", font_size=24, color=ORANGE)
        self.add_fixed_in_frame_mobjects(label_A, label_B)
        label_A.to_edge(LEFT).shift(UP * 2.5)
        label_B.to_edge(RIGHT).shift(UP * 2.5)

        self.play(Create(cube_A), Create(cube_B), FadeIn(label_A), FadeIn(label_B), run_time=1.5)
        self.wait(0.5)

        # ── Vectores evaluados desde cada elemento en su propia cara x=1 (A) / x=0->offset (B) ──
        # Glifo ligero (línea + punta pequeña) en vez de Arrow3D: Arrow3D con
        # su malla cono+cilindro es demasiado costoso de renderizar con el
        # renderer Cairo para decenas de copias (~30s/frame comprobado).
        def make_vectors(offset_shift, vec_cols, color):
            vecs = VGroup()
            for row in data:
                x, y, z = row[0], row[1], row[2]
                vx, vy, vz = row[vec_cols[0]], row[vec_cols[1]], row[vec_cols[2]]
                p0 = (np.array([x, y, z]) + np.array(offset_shift) - np.array([1.0, 0.5, 0.5])) * SCALE
                v = np.array([vx, vy, vz]) * ARROW_SCALE
                if np.linalg.norm(v) < 1e-6:
                    continue
                p1 = p0 + v
                shaft = Line3D(p0, p1, color=color, thickness=0.015)
                tip = Dot3D(p1, radius=0.035, color=color)
                vecs.add(shaft, tip)
            return vecs

        vectors_A = make_vectors((0, 0, 0), (3, 4, 5), BLUE_D)
        vectors_B = make_vectors((2.2, 0, 0), (6, 7, 8), ORANGE)

        self.add(vectors_A, vectors_B)
        self.wait(1)

        # ── Los elementos se acercan hasta compartir la cara ───────────
        self.play(
            cube_B.animate.shift(LEFT * (2.2 - 1.0) * SCALE),
            vectors_B.animate.shift(LEFT * (2.2 - 1.0) * SCALE),
            FadeOut(label_A), FadeOut(label_B),
            run_time=2,
        )
        self.wait(0.5)

        match_text = Text(
            "max |N_A - N_B| en la cara compartida ~ 8e-14  (verificado con BasisEval.jl)",
            font_size=22, color=GREEN,
        ).to_edge(DOWN)
        self.add_fixed_in_frame_mobjects(match_text)
        self.play(Write(match_text))
        self.wait(1)

        self.begin_ambient_camera_rotation(rate=0.15)
        self.wait(2)
