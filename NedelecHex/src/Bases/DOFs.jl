module DOFs

using LinearAlgebra
using ..GaussQuadrature: gauss_segment, gauss_quadrilateral, gauss_hex
using ..CoordinateTransforms: transform_2d_to_3d
using ..ReferenceHex: NODES, EDGES

export evaluate_edge_dof, evaluate_face_dof, evaluate_volume_dof,
    evaluate_all_dofs

# Polinomios de Lagrange en [0,1] para las aristas
const Q_EDGE = (
    t -> 1.0 - t,
    t -> t
)

# Vectores tangentes a cada cara (ec. 7): 6 caras x 4 vectores
const Q_FACE = (
    # Cara 1 (z=0): varía x, y
    (
        pt -> [1.0 - pt[1], 0.0, 0.0],
        pt -> [pt[1], 0.0, 0.0],
        pt -> [0.0, 1.0 - pt[2], 0.0],
        pt -> [0.0, pt[2], 0.0],
    ),
    # Cara 2 (z=1): varía x, y
    (
        pt -> [1.0 - pt[1], 0.0, 0.0],
        pt -> [pt[1], 0.0, 0.0],
        pt -> [0.0, 1.0 - pt[2], 0.0],
        pt -> [0.0, pt[2], 0.0],
    ),
    # Cara 3 (y=0): varía x, z
    (
        pt -> [1.0 - pt[1], 0.0, 0.0],
        pt -> [pt[1], 0.0, 0.0],
        pt -> [0.0, 0.0, 1.0 - pt[3]],
        pt -> [0.0, 0.0, pt[3]],
    ),
    # Cara 4 (x=1): varía y, z
    (
        pt -> [0.0, 1.0 - pt[2], 0.0],
        pt -> [0.0, pt[2], 0.0],
        pt -> [0.0, 0.0, 1.0 - pt[3]],
        pt -> [0.0, 0.0, pt[3]],
    ),
    # Cara 5 (y=1): varía x, z
    (
        pt -> [1.0 - pt[1], 0.0, 0.0],
        pt -> [pt[1], 0.0, 0.0],
        pt -> [0.0, 0.0, 1.0 - pt[3]],
        pt -> [0.0, 0.0, pt[3]],
    ),
    # Cara 6 (x=0): varía y, z
    (
        pt -> [0.0, 1.0 - pt[2], 0.0],
        pt -> [0.0, pt[2], 0.0],
        pt -> [0.0, 0.0, 1.0 - pt[3]],
        pt -> [0.0, 0.0, pt[3]],
    ),
)

# Vectores volumétricos (ec. 8): Lagrange en las tres coordenadas
const Q_VOL = (
    pt -> [1.0 - pt[1], 0.0, 0.0],
    pt -> [pt[1], 0.0, 0.0],
    pt -> [0.0, 1.0 - pt[2], 0.0],
    pt -> [0.0, pt[2], 0.0],
    pt -> [0.0, 0.0, 1.0 - pt[3]],
    pt -> [0.0, 0.0, pt[3]],
)

# g_e(u) = ∫₀¹ (u(γ(t)) · τ̂) · q(t) · |J_e| dt
function evaluate_edge_dof(monomial, edge_idx::Int)
    n1, n2 = EDGES[edge_idx]
    p1 = Float64.(NODES[n1])
    p2 = Float64.(NODES[n2])

    J_edge = p2 - p1
    tau = J_edge / norm(J_edge)
    len = norm(J_edge)

    weights, points = gauss_segment(4)
    # El eltype (Float64 o ComplexF64) se deduce sondeando el campo en p1
    dofs = zeros(typeof(sum(monomial(p1) .* tau)), 2)

    for qi in 1:2
        for k in eachindex(weights)
            t = points[k]
            pt3D = p1 + J_edge * t
            u = monomial(pt3D)
            q = Q_EDGE[qi](t)
            dofs[qi] += weights[k] * sum(u .* tau) * q * len
        end
    end
    return dofs
end

# g_f(u) = ∫_f u · q dS   (dot directo, q ya es tangente a la cara)
function evaluate_face_dof(monomial, face_idx::Int)
    weights, points = gauss_quadrilateral(4)
    det_J2D = 1.0    # hacerlo más general, en vez de Hardcodear 1
    pt0 = transform_2d_to_3d(points[1, :], face_idx)
    dofs = zeros(typeof(sum(monomial(pt0) .* Q_FACE[face_idx][1](pt0))), 4)

    for qi in 1:4
        for k in eachindex(weights)
            pt2D = points[k, :]
            pt3D = transform_2d_to_3d(pt2D, face_idx)
            u = monomial(pt3D)
            q = Q_FACE[face_idx][qi](pt3D)
            dofs[qi] += weights[k] * sum(u .* q) * det_J2D
        end
    end
    return dofs
end

# g_H(u) = ∫_H u · q dV
function evaluate_volume_dof(monomial)
    weights, points = gauss_hex(4)
    det_J = 1.0    # jacobiano = I₃ en el elemento de referencia
    pt0 = points[1, :]
    dofs = zeros(typeof(sum(monomial(pt0) .* Q_VOL[1](pt0))), 6)

    for qi in 1:6
        for k in eachindex(weights)
            pt = points[k, :]
            u = monomial(pt)
            q = Q_VOL[qi](pt)
            dofs[qi] += weights[k] * sum(u .* q) * det_J
        end
    end
    return dofs
end

# Aplica los 54 funcionales a un campo sobre el elemento de referencia, en el
# siguiente orden (12 aristas × 2, 6 caras × 4, 6 de volumen): el mismo orden
# de filas de la matriz G (Fig. 2 del paper).
function evaluate_all_dofs(field)
    edges = reduce(vcat, evaluate_edge_dof(field, e) for e in 1:12)
    faces = reduce(vcat, evaluate_face_dof(field, f) for f in 1:6)
    return vcat(edges, faces, evaluate_volume_dof(field))
end

end