module GlobalDOFs
# Hacer que los elementos que compartan aristas y caras tengan el mismo Dof global 
# (+ @Edge/generateDOF.m, @Face/generateDOF.m, @Element/generateVolumeDOF.m;
# orden local en @Element/getDOF.m).
#También incluye la equivalencia de +assembly/imposeLocalDirichletMatrix.m.

using ..Mesh: HexMesh, edge_polyface_type, face_polyface_type, has_boundary_tags
using ..ReferenceHex: EDGES, FACES

export DOFMap, build_dof_map, n_edge_dofs, n_face_dofs, n_vol_dofs,
    dirichlet_dofs, constrained_dofs, element_dirichlet_mask

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
        # orden local en @Element/getDOF.m:9-18)
        n_dofs += 6
        for i in 1:6
            conn[e, 48+i] = n_dofs - 6 + i
        end
    end

    return DOFMap(n_dofs, conn)
end

# ---------------------------------------------------------------------------
# Condiciones de contorno esenciales (Dirichlet)
# ---------------------------------------------------------------------------
# Una condición esencial fija el valor de ciertos DOFs (los elimina del
# sistema). El código NO está atado a un único tipo: recorre las entidades de
# la malla, consulta el `type` de su PolyFace y despacha la acción según ese
# tipo. Hoy el único tipo esencial es PEC ("PerfectE"), pero el mecanismo es
# extensible: añadir otro tipo esencial (p. ej. un plano de simetría) es
# añadirlo a _is_essential_type y a los helpers de despacho.
#
# No se mira coordenadas: se pregunta a la malla qué condición de contorno
# para PEC parcial, tabiques metálicos interiores, geometrías no alineadas con
# los ejes o problemas con varios tipos de contorno: lo que cambia es el
# etiquetado, no esta función. Las etiquetas las pone el mallador (box_mesh
# vía grupos físicos de Gmsh, load_matlab_mesh vía sus CSV) o
# tag_exterior_boundary!.
#
# antes de ensamblar; aquí se listan los índices globales)

# Tipos de condición de contorno que imponen una restricción esencial.
# Extensible: añade aquí un tipo nuevo y define su acción en los helpers.
_is_essential_type(type::String) = type == "PerfectE"

# DOFs de arista (loc. 1-24) que fija una condición esencial de tipo `type`.
# (PEC: los 2 DOFs de la arista; otro tipo esencial añadiría su propia rama.)
function _essential_edge_dofs!(dir_set::Set{Int}, type::String, dofs, k::Int)
    _is_essential_type(type) || return
    base = 2 * (k - 1)
    push!(dir_set, dofs[base+1])
    push!(dir_set, dofs[base+2])
end

# DOFs de cara (loc. 25-48) que fija una condición esencial de tipo `type`.
# (PEC: los 4 DOFs de la cara; otro tipo esencial añadiría su propia rama.)
function _essential_face_dofs!(dir_set::Set{Int}, type::String, dofs, f::Int)
    _is_essential_type(type) || return
    base = 24 + 4 * (f - 1)
    for i in 1:4
        push!(dir_set, dofs[base+i])
    end
end

# Índices LOCALES (1-54) de los DOFs esenciales del elemento `e`, consultando
# `honor` decide qué tipos de condición se tienen en cuenta. Devuelve un
# Vector{Bool} de longitud 54 (true = ese DOF local es esencial).
function _element_essential_mask(mesh::HexMesh, e::Int, honor::Function)
    gn = mesh.elements[e, :]
    mask = falses(54)

    # - DOF de arista (loc. 1-24): los dos DOFs si la arista es esencial.
    for (k, (ln1, ln2)) in enumerate(EDGES)
        t = edge_polyface_type(mesh, gn[ln1], gn[ln2])
        honor(t) || continue
        base = 2 * (k - 1)
        mask[base+1] = true
        mask[base+2] = true
    end

    # - DOF de cara (loc. 25-48): los cuatro DOFs si la cara es esencial.
    for (f, face_lnodes) in enumerate(FACES)
        gfn = ntuple(i -> gn[face_lnodes[i]], 4)
        t = face_polyface_type(mesh, gfn)
        honor(t) || continue
        base = 24 + 4 * (f - 1)
        for i in 1:4
            mask[base+i] = true
        end
    end
    # - DOFs de volumen (49-54): nunca son esenciales.

    return mask
end

# Núcleo común: recorre la malla y fija los DOFs de las entidades cuyo tipo
# cumple `honor`. `honor` decide qué tipos de condición se tienen en cuenta.
function _collect_essential_dofs(mesh::HexMesh, dof_map::DOFMap, honor::Function)
    has_boundary_tags(mesh) ||
        error("La malla no tiene etiquetas de contorno: no hay forma de saber " *
            "qué DOFs son Dirichlet. Etiquétala (tag_exterior_boundary!, " *
            "load_matlab_mesh o box_mesh) antes de llamar a esta función.")

    dir_set = Set{Int}()
    n_elem = size(mesh.elements, 1)

    for e in 1:n_elem
        dofs = dof_map.connectivity[e, :]
        mask = _element_essential_mask(mesh, e, honor)
        for i in 1:54
            mask[i] && push!(dir_set, dofs[i])
        end
    end

    return sort!(collect(dir_set))
end

# Máscara (54 Bool) de los DOFs LOCALES esenciales del elemento `e`, con el
# predicado por defecto (PEC). La usa el ensamblaje para imponer Dirichlet por
function element_dirichlet_mask(mesh::HexMesh, e::Int)
    return _element_essential_mask(mesh, e, _is_essential_type)
end

# Unificado: DOFs fijados por TODAS las condiciones esenciales de la malla.
# Es el punto de entrada general: recorre las PolyFace y despacha por tipo.
function constrained_dofs(mesh::HexMesh, dof_map::DOFMap)
    return _collect_essential_dofs(mesh, dof_map, _is_essential_type)
end

# Compatibilidad: DOFs fijados por las condiciones esenciales de UN tipo
# concreto (por defecto PEC). Equivale a constrained_dofs cuando solo hay PEC.
function dirichlet_dofs(mesh::HexMesh, dof_map::DOFMap; type::String="PerfectE")
    return _collect_essential_dofs(mesh, dof_map, t -> t == type)
end

end # module GlobalDOFs

#Aclaracíon de diferencia con Matlab:
#numera los DOFs de arista, cara y volumen recorriendo los cada uno en un bucle dentro de la malla
#se garantiza el no repetir Dofs ya que en la malla ya están las entidades solo una vez
#ya que son compartidos por referencia (son objetos handle)
#En Julia recorremos elemento a elemento y comprobamos lo compartido usando diccionarios
#
