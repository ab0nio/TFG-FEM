import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from manim import *
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_matrix(name):
    return np.loadtxt(os.path.join(DATA_DIR, name), delimiter=",")


class BaseDualTablaII(Scene):
    """
    Idea 2: anima el proceso G -> W = G^-1 usando las matrices REALES
    calculadas por NedelecHex.CoefficientGenerator.get_coefficients()
    (exportadas a data/G_matrix.csv y data/W_matrix.csv).
    """

    def construct(self):
        G = load_matrix("G_matrix.csv")
        W = load_matrix("W_matrix.csv")
        n = G.shape[0]
        assert G.shape == (54, 54) and W.shape == (54, 54)

        title = Text("De los DOFs a la Tabla II: G -> W = G⁻¹", font_size=32).to_edge(UP)
        self.play(Write(title))

        # ── Etiquetas de monomios (izq) y funcionales DOF (der) ──────────
        mono_label = VGroup(
            Text("54 monomios", font_size=24),
            MathTex(r"u_1,\,u_2,\,\ldots,\,u_{54}", font_size=28),
        ).arrange(DOWN).to_edge(LEFT).shift(UP * 2)

        dof_label = VGroup(
            Text("54 funcionales DOF", font_size=24),
            MathTex(r"g_1,\,g_2,\,\ldots,\,g_{54}", font_size=28),
        ).arrange(DOWN).to_edge(RIGHT).shift(UP * 2)

        self.play(FadeIn(mono_label), FadeIn(dof_label))
        self.wait(0.5)

        # ── Heatmap de G: 54x54 cuadrados, color = |G_ij| normalizado ────
        cell = 0.09
        gmax = np.abs(G).max()

        def matrix_heatmap(mat, vmax, color_scale):
            grid = VGroup()
            for i in range(n):
                for j in range(n):
                    v = abs(mat[i, j]) / vmax if vmax > 0 else 0.0
                    sq = Square(side_length=cell, stroke_width=0.3, stroke_color=GRAY_D)
                    sq.set_fill(color_scale, opacity=float(min(v, 1.0)))
                    sq.move_to(np.array([j * cell, -i * cell, 0]))
                    grid.add(sq)
            grid.move_to(ORIGIN)
            return grid

        g_grid = matrix_heatmap(G, gmax, BLUE)
        g_caption = Text("G  (g_i aplicado a u_j)", font_size=22).next_to(g_grid, DOWN)

        self.play(FadeOut(mono_label), FadeOut(dof_label))
        self.play(FadeIn(g_grid), Write(g_caption), run_time=2.5)
        self.wait(1)

        # ── Transformación G -> W = G^-1 ──────────────────────────────────
        arrow_txt = MathTex(r"W = G^{-1}", font_size=40)
        self.play(FadeOut(g_caption))
        self.play(g_grid.animate.shift(LEFT * 2.5))
        arrow = Arrow(LEFT, RIGHT, buff=0.2).next_to(g_grid, RIGHT, buff=0.3)
        arrow_txt.next_to(arrow, UP)
        self.play(Create(arrow), Write(arrow_txt))

        wmax = np.abs(W).max()
        w_grid = matrix_heatmap(W, wmax, GREEN)
        w_grid.next_to(arrow, RIGHT, buff=0.3)
        w_caption = Text("W  (coeficientes enteros)", font_size=22).next_to(w_grid, DOWN)

        self.play(FadeIn(w_grid), Write(w_caption), run_time=2.5)
        self.wait(1)

        # ── Zoom sobre un bloque real 5x5 de W: coeficientes enteros ──────
        block = W[0:5, 0:5]
        rows = [
            [f"{val:.0f}" for val in row]
            for row in block
        ]
        matrix_tex = Matrix(rows, h_buff=1.3).scale(0.7)
        zoom_caption = Text("N₁...N₅ sobre a₁...a₅  (bloque real de W, coincide con la Tabla II)",
                             font_size=20)

        self.play(FadeOut(g_grid), FadeOut(w_caption), FadeOut(arrow), FadeOut(arrow_txt),
                   FadeOut(title))
        self.play(w_grid.animate.scale(0.01).move_to(ORIGIN))
        self.remove(w_grid)

        zoom_group = VGroup(matrix_tex, zoom_caption).arrange(DOWN, buff=0.4).move_to(ORIGIN)
        self.play(FadeIn(zoom_group), run_time=1.5)
        self.wait(2)
