module GlobalDOFs

using ..Mesh: HexMesh
using ..ReferenceHex: EDGES, FACES

export DOFMap, build_dof_map, n_edge_dofs, n_face_dofs, n_vol_dofs,
       dirichlet_dofs

const N_EDGE_DOFS_PER_EDGE = 2
const N_FACE_DOFS_PER_FACE = 4
const N_VOL_DOFS_PER_ELEM  = 6

struct DOFMap
    n_dofs::Int
    connectivity::Matrix{Int}   # n_elements × 54
end

Base.show(io::IO, d::DOFMap) =
    print(io, "DOFMap($(d.n_dofs) global DOFs, $(size(d.connectivity,1)) elements)")

n_edge_dofs(d::DOFMap) = count_dof_block(d, 1, 24)
n_face_dofs(d::DOFMap) = count_dof_block(d, 25, 48)
n_vol_dofs(d::DOFMap)  = count_dof_block(d, 49, 54)
count_dof_block(d, a, b) = length(unique(d.connectivity[:, a:b]))

# Build the global DOF connectivity table for a hex mesh.
#
# Edge DOFs (local 1-24): two DOFs per edge.  When a local edge runs opposite
# to the canonical direction (min→max global node), the two DOFs are swapped
# so the global "weight-(1-t)" and "weight-t" DOFs remain consistent.
#
# Face DOFs (local 25-48): four DOFs per face, shared directly between the two
# elements that touch the face.  This relies on Q_FACE using identical test
# vectors for every opposite-face pair (faces 1↔2, 3↔5, 4↔6): valid for
# axis-aligned structured meshes.
#
# Volume DOFs (local 49-54): one set of 6 per element, never shared.
function build_dof_map(mesh::HexMesh)
    n_elem = size(mesh.elements, 1)
    conn   = zeros(Int, n_elem, 54)
    n_dofs = 0

    edge_dofs = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()   # canonical edge → (g1, g2)
    face_dofs = Dict{NTuple{4,Int},  NTuple{4,Int}}()    # sorted face nodes → (g1..g4)

    for e in 1:n_elem
        gn = mesh.elements[e, :]   # 8 global node indices for this element

        # Edge DOFs
        for (k, (ln1, ln2)) in enumerate(EDGES)
            gn1, gn2 = gn[ln1], gn[ln2]
            key = minmax(gn1, gn2)   # canonical orientation: lower → higher

            if !haskey(edge_dofs, key)
                n_dofs += 2
                edge_dofs[key] = (n_dofs - 1, n_dofs)
            end
            g1, g2 = edge_dofs[key]

            base = 2 * (k - 1)
            if gn1 < gn2            # local matches canonical → direct
                conn[e, base + 1] = g1
                conn[e, base + 2] = g2
            else                    # reversed → swap weight-(1-t) ↔ weight-t
                conn[e, base + 1] = g2
                conn[e, base + 2] = g1
            end
        end

        # Face DOFs
        for (f, face_lnodes) in enumerate(FACES)
            gfn = ntuple(i -> gn[face_lnodes[i]], 4)
            key = NTuple{4,Int}(sort!(collect(gfn)))

            if !haskey(face_dofs, key)
                n_dofs += 4
                face_dofs[key] = (n_dofs - 3, n_dofs - 2, n_dofs - 1, n_dofs)
            end
            gd   = face_dofs[key]
            base = 24 + 4 * (f - 1)
            conn[e, base + 1] = gd[1]
            conn[e, base + 2] = gd[2]
            conn[e, base + 3] = gd[3]
            conn[e, base + 4] = gd[4]
        end

        # Volume DOFs
        n_dofs += 6
        for i in 1:6
            conn[e, 48 + i] = n_dofs - 6 + i
        end
    end

    return DOFMap(n_dofs, conn)
end

# Devuelve los índices globales de DOFs que deben anularse por condición PEC
# (E_tan = 0 en toda la frontera exterior del dominio caja).
#
# Criterio geométrico:
# - DOF de arista (loc. 1-24): los dos DOFs se imponen si ambos nodos de la
# arista pertenecen al mismo plano frontera del dominio.
# - DOF de cara (loc. 25-48): los cuatro DOFs se imponen si los cuatro nodos
# de la cara del malla comparten un plano frontera.
# - DOFs de volumen (49-54): nunca son Dirichlet.
function dirichlet_dofs(mesh::HexMesh, dof_map::DOFMap; tol=1e-10)
    coords  = mesh.nodes        # n_nodes × 3
    n_nodes = size(coords, 1)

    xmin, xmax = extrema(coords[:, 1])
    ymin, ymax = extrema(coords[:, 2])
    zmin, zmax = extrema(coords[:, 3])

    # on_plane[n, p] = true si el nodo n está en el plano de frontera p.
    # Planos: 1=x=xmin, 2=x=xmax, 3=y=ymin, 4=y=ymax, 5=z=zmin, 6=z=zmax
    on_plane = falses(n_nodes, 6)
    for n in 1:n_nodes
        on_plane[n, 1] = abs(coords[n, 1] - xmin) < tol
        on_plane[n, 2] = abs(coords[n, 1] - xmax) < tol
        on_plane[n, 3] = abs(coords[n, 2] - ymin) < tol
        on_plane[n, 4] = abs(coords[n, 2] - ymax) < tol
        on_plane[n, 5] = abs(coords[n, 3] - zmin) < tol
        on_plane[n, 6] = abs(coords[n, 3] - zmax) < tol
    end

    dir_set = Set{Int}()
    n_elem  = size(mesh.elements, 1)

    for e in 1:n_elem
        gn   = mesh.elements[e, :]
        dofs = dof_map.connectivity[e, :]

        # Aristas: DOFs locales 1-24 (2 por arista)
        for (k, (ln1, ln2)) in enumerate(EDGES)
            shared = on_plane[gn[ln1], :] .& on_plane[gn[ln2], :]
            if any(shared)
                base = 2 * (k - 1)
                push!(dir_set, dofs[base + 1])
                push!(dir_set, dofs[base + 2])
            end
        end

        # Caras: DOFs locales 25-48 (4 por cara)
        for (f, face_lnodes) in enumerate(FACES)
            shared = reduce(.&, [on_plane[gn[ln], :] for ln in face_lnodes])
            if any(shared)
                base = 24 + 4 * (f - 1)
                for i in 1:4
                    push!(dir_set, dofs[base + i])
                end
            end
        end
        # DOFs de volumen (49-54): nunca Dirichlet
    end

    return sort!(collect(dir_set))
end

end # module GlobalDOFs
