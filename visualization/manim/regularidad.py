import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from manim import *

config.background_color = "#2b2b2b"  # gris oscuro (convención del proyecto)

# ─── Escena: regularidad / formulación fuerte vs débil ───────────────────────
# Ilustra por qué la formulación fuerte (eq:forma-fuerte) exige que la solución
# sea diferenciable dos veces en todo el dominio, mientras que las soluciones
# FEM son polinomios a trozos continuos (C⁰) cuyas derivadas NO están bien
# definidas en las fronteras entre elementos. Eso motiva pasar a la formulación
# débil (integral), que solo pide solución y rotacional de cuadrado integrable.
#
#   → Memoria: subsec:mef-sistema (Métodos), placeholder en memoria_tfg.tex:710.
#   → Exportar como regularidad.png a memoria/figuras/.
#
# 2D y barato (solo Line/Dot/MathTex/Text), según el criterio de IDEAS.md.
#
# Render:  manim -pql regularidad.py RegularidadFormaDebil
# (después mover el PNG a memoria/figuras/regularidad.png)


class RegularidadFormaDebil(Scene):
    """Dos paneles alineados verticalmente:

    Panel superior (u):  solución exacta suave (seno) frente a la aproximación
        FEM poligonal a trozos, con "codos" en las fronteras de elemento.
    Panel inferior (u'): la derivada de la exacta es continua; la derivada de la
        aproximación FEM es constante a trozos y SALTA en cada frontera.
    """

    def construct(self):
        # ── Datos del problema 1D ────────────────────────────────────────────
        # Dominio x ∈ [0, 1], malla uniforme de N elementos.
        N = 6
        xs = np.linspace(0.0, 1.0, N + 1)

        # Solución exacta suave: u(x) = sin(pi * x)  (cumple la forma fuerte).
        def u_exact(x):
            return np.sin(np.pi * x)

        # Aproximación FEM: interpolante lineal a trozos de u_exact en los nodos.
        u_fem_nodes = u_exact(xs)

        # ── Ejes del panel superior (u) ──────────────────────────────────────
        ax_u = Axes(
            x_range=[0, 1, 0.25],
            y_range=[0, 1.1, 0.5],
            x_length=7.0,
            y_length=2.2,
            axis_config={"color": GREY_B, "include_ticks": False},
        ).move_to(UP * 1.6)
        lbl_u = MathTex(r"u(x)", color=WHITE, font_size=28).next_to(
            ax_u, LEFT, buff=0.4)

        # Curva exacta (suave, continua): muestreo fino.
        t_fine = np.linspace(0.0, 1.0, 400)
        exact_curve = ax_u.plot_line_graph(
            x_values=t_fine,
            y_values=u_exact(t_fine),
            line_color=TEAL,
            stroke_width=3.0,
            add_vertex_dots=False,
        )

        # Aproximación FEM: poligonal a trozos (Line por elemento).
        fem_lines = VGroup()
        for i in range(N):
            p0 = ax_u.c2p(xs[i], u_fem_nodes[i])
            p1 = ax_u.c2p(xs[i + 1], u_fem_nodes[i + 1])
            fem_lines.add(Line(p0, p1, color=YELLOW, stroke_width=3.0))

        # Nodos de la malla (donde la derivada salta).
        fem_nodes = VGroup()
        for x, y in zip(xs, u_fem_nodes):
            fem_nodes.add(Dot(ax_u.c2p(x, y), radius=0.05, color=YELLOW))

        # Marcas verticales en las fronteras de elemento (solo interiores).
        boundaries = VGroup()
        for x in xs[1:-1]:
            boundaries.add(
                DashedLine(
                    ax_u.c2p(x, 0), ax_u.c2p(x, 1.1),
                    color=GREY_D, stroke_width=1.5, dash_length=0.06,
                )
            )

        # ── Ejes del panel inferior (u') ─────────────────────────────────────
        ax_du = Axes(
            x_range=[0, 1, 0.25],
            y_range=[-3.6, 3.6, 1.0],
            x_length=7.0,
            y_length=2.2,
            axis_config={"color": GREY_B, "include_ticks": False},
        ).move_to(DOWN * 1.6)
        lbl_du = MathTex(r"u'(x)", color=WHITE, font_size=28).next_to(
            ax_du, LEFT, buff=0.4)

        # Derivada exacta: u'(x) = pi * cos(pi * x)  (continua).
        def du_exact(x):
            return np.pi * np.cos(np.pi * x)

        exact_deriv = ax_du.plot_line_graph(
            x_values=t_fine,
            y_values=du_exact(t_fine),
            line_color=TEAL,
            stroke_width=3.0,
            add_vertex_dots=False,
        )

        # Derivada FEM: constante a trozos = pendiente de cada segmento.
        du_fem_lines = VGroup()
        for i in range(N):
            slope = (u_fem_nodes[i + 1] - u_fem_nodes[i]) / (xs[i + 1] - xs[i])
            p0 = ax_du.c2p(xs[i], slope)
            p1 = ax_du.c2p(xs[i + 1], slope)
            du_fem_lines.add(Line(p0, p1, color=YELLOW, stroke_width=3.0))

        # Marcas de los saltos en la derivada (círculos abiertos en cada extremo
        # del escalón, para indicar que el valor salta en la frontera).
        jump_marks = VGroup()
        for i in range(N):
            slope = (u_fem_nodes[i + 1] - u_fem_nodes[i]) / (xs[i + 1] - xs[i])
            for x in (xs[i], xs[i + 1]):
                jump_marks.add(
                    Circle(radius=0.045, color=YELLOW, stroke_width=2.0).move_to(
                        ax_du.c2p(x, slope))
                )

        # ── Leyenda (con caja de fondo para que se lea bien) ─────────────────
        legend_exact = VGroup(
            Line(ORIGIN, RIGHT * 0.55, color=TEAL, stroke_width=5.0),
            MathTex(r"\text{solución exacta}",
                    color=WHITE, font_size=26),
        ).arrange(RIGHT, buff=0.25)
        legend_fem = VGroup(
            Line(ORIGIN, RIGHT * 0.55, color=YELLOW, stroke_width=5.0),
            MathTex(r"\text{aprox. FEM}",
                    color=WHITE, font_size=26),
        ).arrange(RIGHT, buff=0.25)
        legend_items = VGroup(legend_exact, legend_fem).arrange(
            DOWN, aligned_edge=LEFT, buff=0.22)
        legend_box = SurroundingRectangle(
            legend_items, color=GREY_B, buff=0.25,
            corner_radius=0.1, fill_color=BLACK, fill_opacity=0.55,
            stroke_width=1.5,
        )
        legend = VGroup(legend_box, legend_items).to_corner(UR, buff=0.5)

        # ── Animación ────────────────────────────────────────────────────────
        # Panel superior: primero la solución exacta, luego la FEM encima.
        self.play(Create(ax_u), Write(lbl_u), run_time=1)
        self.play(Create(exact_curve), run_time=1.5)
        self.wait(0.3)
        self.play(FadeIn(boundaries), run_time=0.8)
        self.play(Create(fem_lines), FadeIn(fem_nodes), run_time=1.5)
        self.wait(0.5)

        # Panel inferior: derivadas.
        self.play(Create(ax_du), Write(lbl_du), run_time=1)
        self.play(Create(exact_deriv), run_time=1.5)
        self.wait(0.3)
        self.play(Create(du_fem_lines), FadeIn(jump_marks), run_time=1.5)
        self.wait(0.5)

        # Leyenda y cierre.
        self.play(FadeIn(legend), run_time=1)
        self.wait(2)


if __name__ == "__main__":
    # Renderizado directo:  manim -pql regularidad.py RegularidadFormaDebil
    pass
