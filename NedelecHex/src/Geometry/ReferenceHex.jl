module ReferenceHex
#= 
Traducción de: 
- geometry.Element (constructor 'referenceElement')
- getEdgeByLocalNodes.m
- getFaceByLocalEdges.m
Propósito: Definir la topología fija del hexaedro unidad [0,1]^3.
=#


const NODES = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
]

const EDGES = [
    (1, 2), (2, 3), (4, 3), (1, 4),
    (5, 6), (6, 7), (8, 7), (5, 8),
    (1, 5), (2, 6), (3, 7), (4, 8)
]

const FACE_NORMALS = [
    [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0],
    [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]
]

end # module ReferenceHex
