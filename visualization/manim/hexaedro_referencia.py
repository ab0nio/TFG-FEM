
# Hexaedro referencia fix · PY
"""Hexaedro de referencia H = [0,1]^3 con los 54 DOFs en sus posiciones.
"""
 
import os
import sys
 
sys.path.insert(0, os.path.dirname(__file__))
 
import numpy as np
from manim import *
 
config.background_color = "#2b2b2b"  # gris oscuro
 
REF_SIDE = 4.2  # lado del hexaedro en pantalla (mas grande)
 
# -- Topologia (base 0) -------------------------------------------------------
EDGES = [
    (0, 1), (1, 2), (3, 2), (0, 3),
    (4, 5), (5, 6), (7, 6), (4, 7),
    (0, 4), (1, 5), (2, 6), (3, 7),
]
NODES = np.array([
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
], dtype=float)
 
 
def ref(p):
    """[0,1]^3 -> escena, centrado en el origen."""
    return (np.asarray(p, dtype=float) - 0.5) * REF_SIDE
 
 
def glyph(center, direction, color, length=0.40):
    """Flecha con punta para los DOFs."""
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-12:
        return VGroup()
    d = d / n * length
    p0, p1 = center - d / 2, center + d / 2
    return Arrow(p0, p1, color=color, stroke_width=3.2, buff=0,
                 max_tip_length_to_length_ratio=0.35,
                 tip_shape=StealthTip)


def dot_glyph(center, color, radius=0.09):

    return Dot3D(point=center, radius=radius, color=color)
 
 
class HexaedroReferencia(ThreeDScene):
    def construct(self):
        corners = [ref(n) for n in NODES]
 
        # -- Marco fino alrededor del diagrama --------------------------------
        frame = Rectangle(width=13.0, height=7.0, color=BLACK, stroke_width=1.5)
 
        # -- Wireframe del cubo (lineas mas finas y discretas) -----------------
        wire = VGroup()
        for a, b in EDGES:
            wire.add(Line(corners[a], corners[b], color=GREY_C, stroke_width=3.0))
 
        # -- 24 DOFs de arista (rojo): 2 por arista ----------------------------
        edge_glyphs = VGroup()
        for a, b in EDGES:
            p1, p2 = corners[a], corners[b]
            tau = p2 - p1
            for t in (0.3, 0.7):
                c = p1 + tau * t
                edge_glyphs.add(glyph(c, tau, RED, length=0.40))
 
        # -- 24 DOFs de cara (azul): 4 por cara ---------------------------------
        face_defs = [
            (2, 0.0, 0, 1), (2, 1.0, 0, 1),   # z=0, z=1: varian x,y
            (1, 0.0, 0, 2),                    # y=0: varian x,z
            (0, 1.0, 1, 2),                    # x=1: varian y,z
            (1, 1.0, 0, 2),                    # y=1: varian x,z
            (0, 0.0, 1, 2),                    # x=0: varian y,z
        ]
        face_glyphs = VGroup()
        for fixed_ax, fixed_val, u_ax, v_ax in face_defs:
            for ax, positions in ((u_ax, [(1 / 3, 1 / 2), (2 / 3, 1 / 2)]),
                                  (v_ax, [(1 / 2, 1 / 3), (1 / 2, 2 / 3)])):
                for pu, pv in positions:
                    p = np.full(3, 0.5)
                    p[fixed_ax], p[u_ax], p[v_ax] = fixed_val, pu, pv
                    d = np.zeros(3)
                    d[ax] = 1.0
                    c = ref(p)
                    face_glyphs.add(glyph(c, d, BLUE, length=0.36))
 
        # -- 6 DOFs de volumen (verde) -------------------------------------------
        vol_glyphs = VGroup()
        perp_axes = {0: np.array([0, 1, 0]), 1: np.array([0, 0, 1]), 2: np.array([1, 0, 0])}
        for ax in range(3):
            for t in (1 / 3, 2 / 3):
                p = np.full(3, 0.5)
                p[ax] = t
                # offset perpendicular (un sentido para t=1/3, el otro para
                # t=2/3) para separar visualmente los 6 glyphs.
                sign = -1 if t < 0.5 else 1
                p_offset = ref(p) + sign * 0.55 * perp_axes[ax]
                c = p_offset
                vol_glyphs.add(dot_glyph(c, GREEN))
 
        # -- Ejes coordenados (blanco) -------------------------------------------  
        origin = ref([0.0, 1.0, 0.0])
        axis_len = REF_SIDE * 1.2
        ax_zeta = Arrow(origin, origin + np.array([0, 0, axis_len]),
                        color=WHITE, stroke_width=5, buff=0)
        ax_eta = Arrow(origin, origin + np.array([0, -axis_len, 0]),
                       color=WHITE, stroke_width=5, buff=0)
        ax_xi = Arrow(origin, origin + np.array([axis_len, 0, 0]),
                      color=WHITE, stroke_width=5, buff=0)
        zeta_lbl = MathTex(r"\zeta", color=WHITE).scale(1.2).move_to(
            origin + np.array([0, 0, axis_len + 0.4]))
        eta_lbl = MathTex(r"\xi", color=WHITE).scale(1.2).move_to(
            origin + np.array([0, -axis_len - 0.4, 0]))
        xi_lbl = MathTex(r"\eta", color=WHITE).scale(1.2).move_to(
            origin + np.array([axis_len + 0.4, 0, 0]))

        # Camara 
        self.set_camera_orientation(phi=70 * DEGREES, theta=-55 * DEGREES, zoom=0.72)

        # Desplaza todo el contenido 3D un poco a la derecha para dejar hueco
        # a la leyenda de colores en el margen izquierdo.
        shift_right = np.array([2.2, 0, 0])
        wire.shift(shift_right)
        edge_glyphs.shift(shift_right)
        face_glyphs.shift(shift_right)
        vol_glyphs.shift(shift_right)
        ax_zeta.shift(shift_right)
        ax_eta.shift(shift_right)
        ax_xi.shift(shift_right)
        zeta_lbl.shift(shift_right)
        eta_lbl.shift(shift_right)
        xi_lbl.shift(shift_right)

        # Etiquetas de los ejes: 
        self.add_fixed_orientation_mobjects(zeta_lbl, eta_lbl, xi_lbl)
        self.add_fixed_in_frame_mobjects(frame)

        # -- Leyenda de colores (fija a la camara, margen izquierdo) ------------
        legend = VGroup()
        legend_items = [
            (RED, r"\text{DOFs de arista (24)}", False),
            (BLUE, r"\text{DOFs de cara (24)}", False),
            (GREEN, r"\text{DOFs de volumen (6)}", True),
        ]
        for i, (color, text, is_dot) in enumerate(legend_items):
            y = 2.2 - i * 0.9
            if is_dot:
                marker = Dot(LEFT * 5.6 + UP * y, radius=0.12, color=color)
            else:
                marker = Arrow(LEFT * 5.9 + UP * y, LEFT * 5.3 + UP * y,
                               color=color, stroke_width=3.2, buff=0,
                               max_tip_length_to_length_ratio=0.35,
                               tip_shape=StealthTip)
            label = MathTex(text, color=WHITE).scale(0.7).next_to(marker, RIGHT, buff=0.3)
            legend.add(VGroup(marker, label))
        self.add_fixed_in_frame_mobjects(legend)

        # Orden de anadido: primero ejes, luego el wireframe, luego los glyphs.
        self.add(ax_zeta, ax_eta, ax_xi)
        self.add(wire)
        self.add(edge_glyphs, face_glyphs, vol_glyphs)
 
        self.wait(1)
 
