module ReferenceHex
# Topología fija del hexaedro unidad [0,1]^3.

const NODES = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
]

# Cada arista apunta en el sentido creciente del eje (ξ, η o ζ) al que es
# paralela: no "de índice de nodo menor a mayor".
const EDGES = [
    (1, 2), (2, 3), (4, 3), (1, 4),
    (5, 6), (6, 7), (8, 7), (5, 8),
    (1, 5), (2, 6), (3, 7), (4, 8)
]

const FACE_NORMALS = [
    [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0],
    [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]
]

# Índices de nodo locales (en NODES) para cada una de las 6 caras.
# Cara 1 (z=0), 2 (z=1), 3 (y=0), 4 (x=1), 5 (y=1), 6 (x=0).
const FACES = (
    (1, 2, 3, 4),
    (5, 6, 7, 8),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (4, 3, 7, 8),
    (1, 4, 8, 5),
)

end # module ReferenceHex
