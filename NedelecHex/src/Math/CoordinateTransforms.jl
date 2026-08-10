module CoordinateTransforms
#= 
Traducción de: +math/transform2Dto3DInReferenceElement.m
Propósito: Mapear puntos de cuadratura 2D (u,v) a las caras 3D del hexaedro de referencia.
=#

export transform_2d_to_3d

# MATLAB: +math/transform2Dto3DInReferenceElement.m:10-25
function transform_2d_to_3d(pt2D, face_idx)
    u, v = pt2D
    face_idx == 1 && return [u, v, 0.0]
    face_idx == 2 && return [u, v, 1.0]
    face_idx == 3 && return [u, 0.0, v]
    face_idx == 4 && return [1.0, u, v]
    face_idx == 5 && return [u, 1.0, v]
    face_idx == 6 && return [0.0, u, v]
    error("Face index out of bounds (1-6)")
end
end