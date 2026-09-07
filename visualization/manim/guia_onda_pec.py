import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from manim import *

# ─── Cavidad semillena de dieléctrico rodeada de PEC (imagen estática) ───────
# Cavidad de sección a x b, propagación a lo largo de Z (vertical en la escena
# original, aquí horizontal). Dos medios: vacío (ε0) en la mitad superior
# (y>0) y dieléctrico (εr) en la mitad inferior (y<0). Las paredes laterales
# son conductor eléctrico perfecto (PEC), dibujadas como una capa gruesa de
# color marrón/naranja: n x E = 0.
A, B, L = 1.6, 1.6, 5.5   # semiancho x, semialto y, semilongitud z (escena)
TH = 0.18                 # grosor de la capa PEC

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


def box_corners(x0, x1, y0, y1, z0, z1):
    return [np.array([x0 + (x1 - x0) * n[0],
                      y0 + (y1 - y0) * n[1],
                      z0 + (z1 - z0) * n[2]]) for n in UNIT_NODES]


def wireframe(corners, color, stroke_width=2.5, opacity=1.0):
    return VGroup(*[
        Line(corners[i], corners[j], color=color,
             stroke_width=stroke_width, stroke_opacity=opacity)
        for i, j in EDGES
    ])


def solid_box(corners, color, opacity=0.9, stroke_width=0):
    """Hexaedro relleno (volumen interior de un medio), desplazado al centro
    de las esquinas dadas. Opacidad alta para evitar artefactos de
    transparencia en el renderizado 3D."""
    box = Prism(dimensions=[
        np.linalg.norm(corners[1] - corners[0]),
        np.linalg.norm(corners[3] - corners[0]),
        np.linalg.norm(corners[4] - corners[0]),
    ], fill_color=color, fill_opacity=opacity, stroke_width=stroke_width)
    box.move_to(np.mean([c for c in corners], axis=0))
    return box


def billboard(mob, camera_pos, up=UP, rot=0, tilt=0):
    """Orienta un texto 3D para que mire a la cámara y se lea recto (sin
    inclinación), manteniendo su posición 3D. camera_pos: posición de la cámara.
    rot: giro adicional (en radianes) alrededor del eje de visión.
    tilt: inclinación (en radianes) alrededor del eje horizontal de pantalla
    (positivo = hacia arriba)."""
    center = mob.get_center()
    view = np.array(camera_pos, dtype=float) - center
    view = view / np.linalg.norm(view)
    # "arriba" de pantalla = proyección del eje z sobre el plano ⊥ a view
    su = np.array(up, dtype=float) - np.dot(np.array(up, dtype=float), view) * view
    su = su / np.linalg.norm(su)
    right = np.cross(su, view)
    right = right / np.linalg.norm(right)
    # Matriz que lleva la base local (x,y,z) a (right, su, view)
    R = np.column_stack([right, su, view])
    mob.shift(-center)
    mob.apply_matrix(R)
    if rot:
        mob.rotate(rot, axis=view)
    if tilt:
        mob.rotate(tilt, axis=right)
    mob.shift(center)
    return mob


class GuiaOndaPEC(ThreeDScene):
    """
    Imagen estática de la guía de onda rectangular rodeada de PEC.
    Renderizar con:  manim render -ql -s guia_onda_pec.py GuiaOndaPEC
    (el flag -s guarda el último fotograma como PNG).
    """

    def construct(self):
        # Cámara fija, sin rotaciones: vista de tres cuartos de la guía.
        phi, theta = 75 * DEGREES, -155 * DEGREES
        self.set_camera_orientation(phi=phi, theta=theta)
        # Aleja la cámara (zoom out) para que la guía completa entre en el
        # encuadre sin cortarse.
        self.move_camera(zoom=0.9)
        cam_pos = np.array([np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)]) * 20.0

        # ══ 1. Volumen de la cavidad (dos medios) ═════════════════════════
        title = Text("Cavidad semillena de dieléctrico rodeada de PEC",
                     font_size=34).to_corner(UL)
        self.add_fixed_in_frame_mobjects(title)

        # Losa dieléctrica (εr) en la mitad inferior (z<0); el vacío (ε0)
        # queda en la mitad superior (z>0), sin relleno. La propagación es a
        # lo largo de x (longitud 2L).
        diel = solid_box(box_corners(-L, L, -A, A, -B, 0), BLUE_D, 0.9)
        guide = wireframe(box_corners(-L, L, -A, A, -B, B), WHITE, 3)

        # Etiquetas ε0 (vacío) y εr (dieléctrico) en el plano del puerto de
        # entrada de la onda (cara abierta x = -L, sin PEC), orientadas a la
        # cámara para leerse rectas.
        eps0 = MathTex(r"\varepsilon_0", font_size=62, color=WHITE)
        eps0.move_to([-L, 0, B / 2])
        billboard(eps0, cam_pos, rot=270 * DEGREES)
        epsr = MathTex(r"\varepsilon_r", font_size=62, color=WHITE)
        epsr.move_to([-L, 0, -B / 2])
        billboard(epsr, cam_pos, rot=270 * DEGREES)

        # ══ 2. Capa gruesa de PEC ═════════════════════════════════════════
        PEC_COLOR = ManimColor("#8B5A2B")  # marrón
        wall_defs = [
            ("y = +a", [-L, L, A, A + TH, -B, B]),
            ("y = -a", [-L, L, -A - TH, -A, -B, B]),
            ("z = +b", [-L, L, -A, A, B, B + TH]),
            ("z = -b", [-L, L, -A, A, -B - TH, -B]),
        ]
        walls = VGroup()
        for _name, (x0, x1, y0, y1, z0, z1) in wall_defs:
            c = box_corners(x0, x1, y0, y1, z0, z1)
            wall = Prism(dimensions=[
                np.linalg.norm(c[1] - c[0]),
                np.linalg.norm(c[3] - c[0]),
                np.linalg.norm(c[4] - c[0]),
            ], fill_color=PEC_COLOR, fill_opacity=0.5, stroke_width=0)
            wall.move_to(np.mean([c[i] for i in range(8)], axis=0))
            walls.add(wall)

        # Etiqueta "PEC" sobre la pared lateral derecha (y = +a), en el borde
        # del puerto de entrada de la onda (x = -L), orientada a la cámara
        # para leerse recta.
        pec_label = Text("PEC", font_size=40, color=WHITE)
        pec_label.move_to([L / 2, TH, -L / 3])
        billboard(pec_label, cam_pos, rot=310 * DEGREES)

        # Composición final (todo añadido directamente, sin animaciones).
        self.add(guide, diel, eps0, epsr, walls, pec_label)
        self.wait(1)