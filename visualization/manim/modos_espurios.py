import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from manim import *

# ─── Escena 8: Modos espurios / núcleo de K ─────────────────────────────────
# Ilustra por qué la matriz de rigidez K tiene un núcleo enorme (campos
# gradiente ⇒ ∇×∇φ = 0), qué le hace eso al espectro (cluster de ~0 + modos
# físicos) y por qué el benchmark usa el solver denso en lugar de shift-invert
# de ARPACK. 2D y barato (solo Line/Dot/MathTex), según el criterio de IDEAS.md.
#
#   → Memoria: sec:nucleo (Métodos) y sec:espurios (Estado del arte).
#   → Exportar como nucleo_K.png a memoria/figuras/.


class NucleoKModosEspurios(Scene):
    """
    Tres actos:
      1. El núcleo: un campo gradiente ∇φ cumple ∇×∇φ = 0 ⇒ K x = 0.
      2. El espectro: cluster de autovalores ~0 (espurios) + modos físicos ~1e5,
         separados por un umbral.
      3. Por qué importa: shift-invert de ARPACK sufre con el cluster; el
         benchmark usa el solver denso.
    """

    def construct(self):
        # ══ 1. El núcleo: ∇×∇φ = 0 ════════════════════════════════════════
        title = Text("Núcleo de la matriz de rigidez  K", font_size=34).to_edge(UP)
        self.play(Write(title))
        self.wait(0.3)

        # Campo gradiente ∇φ (flechas radiales) sobre una rejilla 2D.
        # φ = x² + y²  ⇒  ∇φ = (2x, 2y): campo radial que apunta hacia fuera.
        grid = VGroup()
        arrows = VGroup()
        xs = np.linspace(-2.2, 2.2, 7)
        ys = np.linspace(-1.6, 1.6, 5)
        for x in xs:
            for y in ys:
                grid.add(Dot([x, y, 0], radius=0.03, color=GRAY_D))
                g = np.array([2 * x, 2 * y, 0.0])
                n = np.linalg.norm(g)
                if n < 1e-9:
                    continue
                g = g / n * 0.35
                arrows.add(
                    Line([x, y, 0], [x + g[0], y + g[1], 0],
                         color=YELLOW, stroke_width=3)
                )
        grid.move_to(LEFT * 3.2)
        arrows.move_to(LEFT * 3.2)

        lbl_grad = MathTex(r"\mathbf{u} = \nabla\phi", font_size=36, color=YELLOW)
        lbl_grad.next_to(grid, DOWN, buff=0.5)
        self.play(FadeIn(grid), Create(arrows, lag_ratio=0.05), run_time=2)
        self.play(Write(lbl_grad))
        self.wait(0.5)

        # ∇×∇φ = 0  ⇒  K x = 0
        eq_curl = MathTex(
            r"\nabla\times(\nabla\phi) = \mathbf{0}",
            r"\quad\Rightarrow\quad",
            r"\mathbf{K}\,\mathbf{x} = \mathbf{0}",
            font_size=34,
        )
        eq_curl[0].set_color(YELLOW)
        eq_curl[2].set_color(RED)
        eq_curl.next_to(lbl_grad, DOWN, buff=0.5)
        self.play(Write(eq_curl))
        self.wait(0.8)

        note = Text("Todo campo gradiente está en el núcleo de K",
                    font_size=24, color=GRAY_B)
        note.next_to(eq_curl, DOWN, buff=0.4)
        self.play(FadeIn(note))
        self.wait(1.2)

        # ══ 2. El espectro ════════════════════════════════════════════════
        self.play(
            FadeOut(grid), FadeOut(arrows), FadeOut(lbl_grad),
            FadeOut(eq_curl), FadeOut(note),
            title.animate.set_opacity(0.25),
        )

        spec_title = Text("Espectro de K  (autovalores)", font_size=30).to_edge(UP)
        self.play(Write(spec_title))

        # Eje logarítmico de 1e-12 a 1e6 (8 décadas).
        ax = Axes(
            x_range=[0, 1, 1],
            y_range=[-12, 6, 2],
            x_length=0.4,
            y_length=5.6,
            axis_config={"include_numbers": False, "stroke_width": 2},
        ).to_edge(LEFT, buff=1.2)
        ax.shift(RIGHT * 0.6)
        self.play(Create(ax), run_time=1)

        # Etiquetas de las décadas (potencias de 10).
        for k in range(-12, 7, 2):
            lab = MathTex(r"10^{%d}" % k, font_size=18)
            lab.next_to(ax.c2p(0, k), LEFT, buff=0.15)
            lab.set_opacity(0)
            self.add(lab)
        self.wait(0.2)

        # Cluster espurio: ~10^-9 (del orden de precisión de máquina).
        # 26 modos nulos en la malla 1x1x1 (tabla sec:res-nucleo).
        spurious = VGroup()
        n_spurious = 26
        for i in range(n_spurious):
            y = -9.0 + 0.02 * np.sin(i * 1.7)
            x = 0.5 + 0.35 * np.cos(i * 2.3)
            spurious.add(Dot(ax.c2p(x, y), radius=0.045, color=RED))
        lbl_spur = Text("espurios  ~ 10⁻⁹", font_size=22, color=RED)
        lbl_spur.next_to(ax.c2p(0.5, -9.0), RIGHT, buff=0.3)

        self.play(Create(spurious, lag_ratio=0.02), run_time=2)
        self.play(FadeIn(lbl_spur))
        self.wait(0.5)

        # Modos físicos: ~10^5 (k_c² del modo TE101, del orden de 1e5).
        physical = VGroup()
        for i, y in enumerate([4.9, 5.1, 5.3, 5.6, 5.9]):
            x = 0.5 + 0.3 * np.cos(i * 2.1)
            physical.add(Dot(ax.c2p(x, y), radius=0.05, color=TEAL))
        lbl_phys = Text("físicos  ~ 10⁵", font_size=22, color=TEAL)
        lbl_phys.next_to(ax.c2p(0.5, 5.3), RIGHT, buff=0.3)

        self.play(Create(physical, lag_ratio=0.1), run_time=1.5)
        self.play(FadeIn(lbl_phys))
        self.wait(0.5)

        # Umbral de filtrado: separa espurios de físicos.
        thresh = DashedLine(
            ax.c2p(0.0, -2.0), ax.c2p(1.0, -2.0), color=WHITE, stroke_width=2
        )
        lbl_thresh = Text("umbral", font_size=20, color=WHITE)
        lbl_thresh.next_to(ax.c2p(0.5, -2.0), RIGHT, buff=0.3)
        self.play(Create(thresh), FadeIn(lbl_thresh))
        self.wait(0.6)

        note2 = Text("Separarlos es inmediato: 10⁵ vs 10⁻⁹",
                     font_size=24, color=GRAY_B).to_edge(DOWN)
        self.play(FadeIn(note2))
        self.wait(1.2)

        # ══ 3. Por qué importa: ARPACK shift-invert vs solver denso ═══════
        self.play(
            FadeOut(spurious), FadeOut(physical), FadeOut(lbl_spur),
            FadeOut(lbl_phys), FadeOut(thresh), FadeOut(lbl_thresh),
            FadeOut(note2), FadeOut(spec_title), FadeOut(ax),
            title.animate.set_opacity(0),
        )

        why_title = Text("Consecuencia numérica", font_size=30).to_edge(UP)
        self.play(Write(why_title))

        # El cluster de ~0 hace que shift-invert de ARPACK converja mal.
        eq_arpack = MathTex(
            r"\text{ARPACK shift-invert: }"
            r"(K - \sigma M)^{-1}",
            font_size=30,
        )
        eq_arpack.next_to(why_title, DOWN, buff=0.6)
        self.play(Write(eq_arpack))

        prob = Text(
            "el cluster de ~0 satura la convergencia hacia los modos bajos",
            font_size=22, color=RED,
        ).next_to(eq_arpack, DOWN, buff=0.4)
        self.play(FadeIn(prob))
        self.wait(1)

        # El benchmark usa el solver denso (eigvals completos) para mallas
        # pequeñas/medianas.
        eq_dense = MathTex(
            r"\text{benchmark: solver denso }"
            r"\Rightarrow \text{espectro completo, sin ambigüedad}",
            font_size=28, color=TEAL,
        )
        eq_dense.next_to(prob, DOWN, buff=0.6)
        self.play(Write(eq_dense))
        self.wait(1.5)

        # Cierre: el núcleo es estructura del operador, no un fallo.
        close = Text(
            "El núcleo no es un fallo: es la estructura del rotacional",
            font_size=26, color=GRAY_B,
        ).to_edge(DOWN)
        self.play(FadeIn(close))
        self.wait(2)


if __name__ == "__main__":
    # Renderizado directo:  manim -pql modos_espurios.py NucleoKModosEspurios
    pass
