module GlobalDOFs
# Hacer que los elementos que compartan aristas y caras tengan el mismo Dof global 
# Equivalente en MATLAB: numeración de DOFs en +geometry/@Mesh/generateDOF.m
# (+ @Edge/generateDOF.m, @Face/generateDOF.m, @Element/generateVolumeDOF.m;
# orden local en @Element/getDOF.m).
#También incluye la equivalencia de +assembly/imposeLocalDirichletMatrix.m.

using ..Mesh: HexMesh, edge_polyface_type, face_polyface_type, has_boundary_tags
using ..ReferenceHex: EDGES, FACES

export DOFMap, build_dof_map, n_edge_dofs, n_face_dofs, n_vol_dofs,
    dirichlet_dofs

const N_EDGE_DOFS_PER_EDGE = 2
const N_FACE_DOFS_PER_FACE = 4
const N_VOL_DOFS_PER_ELEM = 6


#Struct de salida con el número de Dofs globales y la tabla de conectividad
#ejemplo de tabla de conectividad n_elements x 54:
# connectivity[e,i] = índice global del DOF local i del elemento e = k 
struct DOFMap
    n_dofs::Int
    connectivity::Matrix{Int}   # n_elementos × 54
end

Base.show(io::IO, d::DOFMap) =
    print(io, "DOFMap($(d.n_dofs) global DOFs, $(size(d.connectivity,1)) elements)")

n_edge_dofs(d::DOFMap) = count_dof_block(d, 1, 24)
n_face_dofs(d::DOFMap) = count_dof_block(d, 25, 48)
n_vol_dofs(d::DOFMap) = count_dof_block(d, 49, 54)
count_dof_block(d, a, b) = length(unique(d.connectivity[:, a:b]))

# Construye la tabla de conectividad global de DOFs para una malla de hexaedros.
function build_dof_map(mesh::HexMesh)
    n_elem = size(mesh.elements, 1)
    conn = zeros(Int, n_elem, 54)
    n_dofs = 0

    edge_dofs = Dict{Tuple{Int,Int},Tuple{Int,Int}}()   # arista canónica → (g1, g2)
    face_dofs = Dict{NTuple{4,Int},NTuple{4,Int}}()    # nodos de cara ordenados → (g1..g4)

    for e in 1:n_elem
        gn = mesh.elements[e, :]   # los 8 índices de nodo globales de este elemento

        # DOFs de arista (locales 1-24): dos DOFs por arista. Cuando una arista local se recorre 
        # en sentido contrario a la dirección canónica (nodo global min→max), los dos DOFs
        # se intercambian para que los pesos t y 1-t sean iguales en todas las aristas
        for (k, (ln1, ln2)) in enumerate(EDGES)
            gn1, gn2 = gn[ln1], gn[ln2]
            key = minmax(gn1, gn2)   # orientación : menor → mayor pq devuelve (min,max)

            if !haskey(edge_dofs, key)#devuelve true si la clave NO está dentro del diccionario
                n_dofs += 2
                edge_dofs[key] = (n_dofs - 1, n_dofs)
            end
            g1, g2 = edge_dofs[key]

            base = 2 * (k - 1)
            if gn1 < gn2            # el local coincide con el canónico → directo
                conn[e, base+1] = g1
                conn[e, base+2] = g2
            else                    # invertido → intercambia peso-(1-t) ↔ peso-t
                conn[e, base+1] = g2
                conn[e, base+2] = g1
            end
        end

        # DOFs de cara (locales 25-48): cuatro DOFs por cara, compartidos directamente
        # entre los dos elementos que tocan la cara. Esto se apoya en que Q_FACE usa
        #vectores de test idénticos para cada par de caras opuestas (caras 1↔2, 3↔5,
        # 4↔6): es válido porque son mallas estructuradas alineadas con los ejes.
        for (f, face_lnodes) in enumerate(FACES)
            gfn = ntuple(i -> gn[face_lnodes[i]], 4)
            key = NTuple{4,Int}(sort!(collect(gfn)))

            if !haskey(face_dofs, key)
                n_dofs += 4
                face_dofs[key] = (n_dofs - 3, n_dofs - 2, n_dofs - 1, n_dofs)
            end
            gd = face_dofs[key]
            base = 24 + 4 * (f - 1)
            conn[e, base+1] = gd[1]
            conn[e, base+2] = gd[2]
            conn[e, base+3] = gd[3]
            conn[e, base+4] = gd[4]
        end

        # DOFs de volumen (locales 49-54): siempre son 6 por elemento.
        # MATLAB: +geometry/@Mesh/generateDOF.m:9-19 (numera aristas/caras/volumen;
        # orden local en @Element/getDOF.m:9-18)
        n_dofs += 6
        for i in 1:6
            conn[e, 48+i] = n_dofs - 6 + i
        end
    end

    return DOFMap(n_dofs, conn)
end

# Devuelve los índices globales de DOFs que deben anularse por condición PEC
# (E_tan = 0 sobre las entidades etiquetadas con `type`).
#
# No mira coordenadas: pregunta a la malla qué condición de contorno tiene cada
# arista y cada cara, igual que MATLAB. Así el mismo código sirve para PEC
# parcial, tabiques metálicos interiores, geometrías no alineadas con los ejes o
# problemas con varios tipos de contorno: lo que cambia es el etiquetado, no
# esta función. Las etiquetas las pone el mallador (box_mesh vía grupos físicos
# de Gmsh, load_matlab_mesh vía exportTags.m) o tag_exterior_boundary!.
#
# MATLAB: +assembly/imposeLocalDirichletMatrix.m (allí se aplica por elemento
# antes de ensamblar; aquí se listan los índices globales)
function dirichlet_dofs(mesh::HexMesh, dof_map::DOFMap; type::String="PerfectE")
    has_boundary_tags(mesh) ||
        error("La malla no tiene etiquetas de contorno: no hay forma de saber " *
              "qué DOFs son Dirichlet. Etiquétala (tag_exterior_boundary!, " *
              "load_matlab_mesh o box_mesh) antes de llamar a dirichlet_dofs.")

    dir_set = Set{Int}()
    n_elem = size(mesh.elements, 1)

    for e in 1:n_elem
        gn = mesh.elements[e, :]
        dofs = dof_map.connectivity[e, :]

        # - DOF de arista (loc. 1-24): los dos DOFs se imponen si la arista
        # pertenece a una polyface del tipo pedido.
        for (k, (ln1, ln2)) in enumerate(EDGES)
            if edge_polyface_type(mesh, gn[ln1], gn[ln2]) == type
                base = 2 * (k - 1)
                push!(dir_set, dofs[base+1])
                push!(dir_set, dofs[base+2])
            end
        end

        # - DOF de cara (loc. 25-48): los cuatro DOFs se imponen si la cara
        # pertenece a una polyface del tipo pedido.
        for (f, face_lnodes) in enumerate(FACES)
            gfn = ntuple(i -> gn[face_lnodes[i]], 4)
            if face_polyface_type(mesh, gfn) == type
                base = 24 + 4 * (f - 1)
                for i in 1:4
                    push!(dir_set, dofs[base+i])
                end
            end
        end
        # - DOFs de volumen (49-54): nunca son Dirichlet.
    end

    return sort!(collect(dir_set))
end

end # module GlobalDOFs

#Aclaracíon de diferencia con Matlab:
#numera los DOFs de arista, cara y volumen recorriendo los cada uno en un bucle dentro de la malla
#se garantiza el no repetir Dofs ya que en la malla ya están las entidades solo una vez
#ya que son compartidos por referencia (son objetos handle)
#En Julia recorremos elemento a elemento y comprobamos lo compartido usando diccionarios
#
