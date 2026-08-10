module Assembly
# Equivalente en MATLAB: +assembly/getMassMatrix.m, +assembly/getStiffMatrix.m,
# +assembly/getGlobalMatrices.m (ensamblaje disperso por formato COO).

using LinearAlgebra
using SparseArrays
using ..GaussQuadrature: gauss_hex
using ..BasisEval: eval_bases_ref, eval_curl_bases_ref,
    compute_jacobian, piola_transform, curl_piola_transform
using ..Mesh: HexMesh, element_nodes, element_epsilon, element_inv_mu
using ..GlobalDOFs: DOFMap

export mass_matrix_element, stiff_matrix_element,
    assemble_global_matrices, apply_dirichlet!, free_dofs,
    load_vector_element, assemble_load_vector

# 2p = 4 → 3 puntos/dirección → 27 puntos de Gauss.
const INTEGRATION_ORDER = 4

# M_ij = ∫ ε N_i·N_j dΩ  (ec. 20 del paper)
# MATLAB: +assembly/getMassMatrix.m:23-32
function mass_matrix_element(W, monomials, nodes; epsilon_r=1.0)
    gaussW, gaussPts = gauss_hex(INTEGRATION_ORDER)
    M = zeros(54, 54)
    for g in eachindex(gaussW)
        pt = @view gaussPts[g, :]
        J = compute_jacobian(nodes, pt)
        N = piola_transform(eval_bases_ref(W, monomials, pt), J)
        M .+= (gaussW[g] * epsilon_r * abs(det(J))) .* (N * N')
    end
    return M
end

# K_ij = ∫ μ⁻¹ curl(N_i)·curl(N_j) dΩ  (ec. 21 del paper)
# MATLAB: +assembly/getStiffMatrix.m:23-31
function stiff_matrix_element(W, curl_monomials, nodes; inv_mu_r=1.0)
    gaussW, gaussPts = gauss_hex(INTEGRATION_ORDER)
    K = zeros(54, 54)
    for g in eachindex(gaussW)
        pt = @view gaussPts[g, :]
        J = compute_jacobian(nodes, pt)
        cN = curl_piola_transform(eval_curl_bases_ref(W, curl_monomials, pt), J)
        K .+= (gaussW[g] * inv_mu_r * abs(det(J))) .* (cN * cN')
    end
    return K
end

# Vector de carga elemental Fᵢ = ∫ Nᵢ·O dΩ, con O el término de forzado fuerte
# del método de las soluciones manufacturadas (O ya evaluado en coordenadas
# físicas). Válido para hexaedros con caras alineadas a los ejes (J diagonal
# y constante en el elemento), que es la única geometría producida por
# box_mesh: la posición física se obtiene por el mapeo afín x = x₀ + J·ξ en
# vez del mapeo isoparamétrico trilineal general.
# Sin equivalente en MATLAB: allí no hay MMS (el único forzado es el problema
# de autovalores de la cavidad, que es homogéneo).
function load_vector_element(W, monomials, O_func, nodes)
    gaussW, gaussPts = gauss_hex(INTEGRATION_ORDER)
    J = compute_jacobian(nodes, [0.5, 0.5, 0.5])
    scale = [J[1, 1], J[2, 2], J[3, 3]]
    x0 = nodes[1]
    F = zeros(54)
    for g in eachindex(gaussW)
        pt = @view gaussPts[g, :]
        x_phys = x0 .+ scale .* pt
        N = piola_transform(eval_bases_ref(W, monomials, pt), J)
        F .+= (gaussW[g] * abs(det(J))) .* (N * O_func(x_phys))
    end
    return F
end

# Ensamblaje global del vector de carga: mismo patrón que
# assemble_global_matrices (suma en los DOFs compartidos), pero sobre un vector.
function assemble_load_vector(mesh::HexMesh, dof_map::DOFMap, W, mono, O_func)
    n_elem = size(mesh.elements, 1)
    F = zeros(dof_map.n_dofs)
    for e in 1:n_elem
        nodes = element_nodes(mesh, e)
        dofs = dof_map.connectivity[e, :]
        Fe = load_vector_element(W, mono, O_func, nodes)
        for i in 1:54
            F[dofs[i]] += Fe[i]
        end
    end
    return F
end

# Ensambla las matrices globales K y M recorriendo todos los elementos.
# Por defecto el material sale de las etiquetas de la malla (mesh.solids +
# mesh.element_solid), como en MATLAB: getGlobalMatrices.m lee
# solids(elements(e).solidId).materialObj. epsilon_r/inv_mu_r se pueden pasar a
# mano (escalar homogéneo o vector de longitud n_elem) para sobreescribirlas.
# MATLAB: +assembly/getGlobalMatrices.m:20-58
function assemble_global_matrices(
    mesh::HexMesh, dof_map::DOFMap, W, mono, curl_mono;
    epsilon_r::Union{Nothing,Float64,AbstractVector{Float64}}=nothing,
    inv_mu_r::Union{Nothing,Float64,AbstractVector{Float64}}=nothing)

    n_dofs = dof_map.n_dofs
    n_elem = size(mesh.elements, 1)
    eps_vec = epsilon_r === nothing ? element_epsilon(mesh) :
              epsilon_r isa Float64 ? fill(epsilon_r, n_elem) : epsilon_r
    mu_vec = inv_mu_r === nothing ? element_inv_mu(mesh) :
             inv_mu_r isa Float64 ? fill(inv_mu_r, n_elem) : inv_mu_r

    # Pre-alocar array COO (rows, cols, value) (cada elemento tiene 54×54 entradas).
    nnz_e = 54 * 54
    total_nz = n_elem * nnz_e
    rows = Vector{Int}(undef, total_nz)
    cols = Vector{Int}(undef, total_nz)
    Mvals = Vector{Float64}(undef, total_nz)
    Kvals = Vector{Float64}(undef, total_nz)

    #bucle de ensamblaje
    for e in 1:n_elem
        nodes = element_nodes(mesh, e)
        dofs = dof_map.connectivity[e, :]

        Me = mass_matrix_element(W, mono, nodes; epsilon_r=eps_vec[e])
        Ke = stiff_matrix_element(W, curl_mono, nodes; inv_mu_r=mu_vec[e])

        base = (e - 1) * nnz_e #offset donde empieza el bloque de este elemento (elemento 1(1-2916) elemento 2(2917-5832) ...)
        for j in 1:54, i in 1:54#recorre las filas i y luego las columnas j del bloque de 54x54 de este elemento
            p = base + (j - 1) * 54 + i
            rows[p] = dofs[i] #traduccción del dof local(i,j) a dof global(dof[i], dof[j])
            cols[p] = dofs[j]
            Mvals[p] = Me[i, j]
            Kvals[p] = Ke[i, j]
        end
    end

    M = sparse(rows, cols, Mvals, n_dofs, n_dofs)#tranasforma COO a CSC que suma los duplicados que apunta a la misma posición (dofs[i], dofs[j]) de la matriz global
    K = sparse(rows, cols, Kvals, n_dofs, n_dofs)

    # Descarta las entradas que son demasiado pequeñas
    droptol!(M, eps() * maximum(abs, nonzeros(M); init=0.0))
    droptol!(K, eps() * maximum(abs, nonzeros(K); init=0.0))

    # droptol! encoge el length de los arrays pero Julia no devuelve la capacidad
    # ya reservada, así que la memoria seguiría ocupada. Reconstruir sí la libera:
    M = sparse(findnz(M)..., n_dofs, n_dofs)
    K = sparse(findnz(K)..., n_dofs, n_dofs)
    return K, M
end

# Devuelve los índices de DOFs libres (complemento de dir_dofs).
free_dofs(n_dofs::Int, dir_dofs::AbstractVector{Int}) = setdiff(1:n_dofs, dir_dofs)

# Impone condiciones Dirichlet: anula filas/columnas y pone 1 en diagonal.
# MATLAB: +assembly/imposeLocalDirichletMatrix.m:46-52 (MATLAB lo hace por
# elemento antes de ensamblar; aquí se hace globalmente después)
function apply_dirichlet!(K::SparseMatrixCSC, M::SparseMatrixCSC,
    dir_dofs::AbstractVector{Int})
    dir_set = Set(dir_dofs) # Set hace que busque en un hash en vez de un vector
    n = size(K, 1)

    for mat in (K, M) #por cada matriz
        for j in 1:n #recorre cada columna
            for p in nzrange(mat, j) # por cada rango de indices no nulos de la columna j
                i = mat.rowval[p] #Busca en memoria la fila i que pertenece a ese elemento en la posición p
                if i in dir_set || j in dir_set #si la fila o columna pertenece a un dof de dirichlet
                    mat.nzval[p] = (i == j && i in dir_set) ? 1.0 : 0.0
                end
            end
        end
        dropzeros!(mat)#elimina los ceros y recomprime a CSC
    end
end

end # module Assembly
