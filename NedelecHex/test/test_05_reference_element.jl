# Prueba BasisEval.jl + Assembly.jl en el elemento de referencia. MATLAB:
# computeJacobianMatrix.m (Jacobiano), getCoefficients.m (base dual),
# getStiffMatrix.m/getMassMatrix.m (Kᵉ/Mᵉ).
@testset "BasisEval.jl: Jacobiano del elemento de referencia" begin
    nodes = NedelecHex.ReferenceHex.NODES

    # Para el hexaedro de referencia [0,1]³ con nodos exactamente en NODES,
    # el mapeo isoparamétrico es la identidad, así que J debe ser I₃ en
    # cualquier punto de referencia.
    explica("El Jacobiano en el Hexaedro de Referencia (al no estar deformado debe ser la identidad)
    J = I₃, det J = 1 y J⁻ᵀ = I₃ en cualquier punto (se muestrean 4, incluidas dos esquinas).")
    for pt in ([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.25, 0.75, 0.5])
        J = NedelecHex.BasisEval.compute_jacobian(nodes, pt)
        @test norm(J - I) < 1e-12
        @test isapprox(det(J), 1.0, atol=1e-12)
        @test norm(inv(J)' - I) < 1e-12
    end

    # Para un hexaedro escalado [0,a]×[0,b]×[0,c], J debe ser diag(a,b,c),
    # det(J) = a·b·c y J⁻ᵀ = diag(1/a,1/b,1/c) (J diagonal ⇒ J⁻ᵀ = J⁻¹).
    explica("Ladrillo [0,2]×[0,3]×[0,0.5]: J = diag(a,b,c), det J = abc = 3.0, J⁻ᵀ = diag(1/a,1/b,1/c)..")
    a, b, c = 2.0, 3.0, 0.5
    scaled_nodes = [n .* [a, b, c] for n in NedelecHex.ReferenceHex.NODES]
    for pt in ([0.5, 0.5, 0.5], [0.1, 0.9, 0.3])
        J_scaled = NedelecHex.BasisEval.compute_jacobian(scaled_nodes, pt)
        @test norm(J_scaled - Diagonal([a, b, c])) < 1e-12
        @test isapprox(det(J_scaled), a * b * c, atol=1e-10)
        @test norm(inv(J_scaled)' - Diagonal([1/a, 1/b, 1/c])) < 1e-10
    end
end
@testset "BasisEval.jl: Propiedad dual de la base y curl" begin
    W, mono, curl_mono = W_REF, MONO, CURL_MONO

    # Propiedad dual: gᵢ(Nⱼ) = δᵢⱼ
    # A diferencia de test_G_matrix.jl (que sólo comprueba G·W=I algebraicamente),
    # aquí se evalúan las funciones de base reales vía eval_bases_ref y se les
    # aplican los funcionales DOF de DOFs.jl directamente, ejercitando el
    # camino completo BasisEval.jl → DOFs.jl de punta a punta.
    D = zeros(54, 54)
    for j in 1:54
        Nj = pt -> NedelecHex.BasisEval.eval_bases_ref(W, mono, pt)[j, :]
        D[:, j] = NedelecHex.DOFs.evaluate_all_dofs(Nj)
    end
    explica("Compruebo que, al aplicar los 54 DOFs a las 54 bases evaluadas en el hexaedro de referencia,
    sale la identidad (propiedad dual gᵢ(Nⱼ) = δᵢⱼ).:
    la matriz D con gᵢ(Nⱼ) se desvía de I₅₄ como máximo $(maximum(abs.(D - I(54)))). .")
    @test isapprox(D, I(54), atol=1e-8)

    # Curl de las funciones de base asociadas al volumen (49-54)
    # En el centro del elemento [0.5,0.5,0.5] el curl se anula por simetría,
    # así que se evalúa en un punto asimétrico.
    explica("Las 6 bases de volumen deben tener curl ≠ 0. (salvo en el centro)")
    pt_asym = [0.2, 0.7, 0.9]
    curlN = NedelecHex.BasisEval.eval_curl_bases_ref(W, curl_mono, pt_asym)
    for j in 49:54
        @test norm(curlN[j, :]) > 1e-6
    end

    # Curl de los modos irrotacionales (espurios)
    # El espacio nulo de Kᵉ = ∫ curlNᵢ·curlNⱼ corresponde a combinaciones
    # "gradiente" de las funciones de base cuyo curl es idénticamente nulo.
    # Se extrae un vector propio del autovalor más pequeño de Kᵉ y se
    # reconstruye el campo curl en varios puntos: debe ser ≈0 en todos ellos.
    nodes = NedelecHex.ReferenceHex.NODES
    Ke = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes)
    eg = eigen(Symmetric(Ke))
    explica("Modos espurios: el autovalor mínimo de Kᵉ es $(eg.values[1]) ≈ 0.
    Demuestra que hay modos con energía de rotacional cero , luego de ese autovalor saca el autovector
    y calcula su curl en 3 puntos para demostrar que es cero y que los modos espurios que se cuelan en la cavidad
    tienen que estar allí")
    @test eg.values[1] < 1e-8         # autovalor nulo (espurio)
    x_null = eg.vectors[:, 1]

    for pt in ([0.5, 0.5, 0.5], [0.2, 0.7, 0.9], [0.6, 0.1, 0.4])
        curl_field = NedelecHex.BasisEval.eval_curl_bases_ref(W, curl_mono, pt)' * x_null
        @test norm(curl_field) < 1e-8
    end
end
@testset "Assembly.jl: Matrices locales Kᵉ, Mᵉ en el elemento de referencia" begin
    W, mono, curl_mono = W_REF, MONO, CURL_MONO
    nodes = NedelecHex.ReferenceHex.NODES

    Ke = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes)
    Me = NedelecHex.Assembly.mass_matrix_element(W, mono, nodes)

    @test size(Ke) == (54, 54)
    @test size(Me) == (54, 54)

    # Simetría
    explica("Kᵉ y Mᵉ simétricas: las formas bilineales c₁, c₂ (ecs. 20-21 del paper) lo son por construcción.")
    @test norm(Ke - Ke') < 1e-12
    @test norm(Me - Me') < 1e-12

    # Mᵉ definida positiva
    λM = eigvals(Symmetric(Me))
    explica("Mᵉ definida positiva (es el producto interno L²): λ_min = $(minimum(λM)) > 0. Kᵉ sólo semidefinida: la energía ∫|∇×u|² ≥ 0 con igualdad en los modos gradiente.")
    @test minimum(λM) > 0

    # Kᵉ semidefinida positiva
    λK = eigvals(Symmetric(Ke))
    @test minimum(λK) > -1e-8

    # Número de condicionamiento (Tabla III del paper: 19 y 35)
    # El paper reporta el condicionamiento de las matrices *precondicionadas*
    # con un escalado diagonal de Jacobi (Mᵖ = D⁻¹ᐟ²MD⁻¹ᐟ², análogo para Kᵖ).
    # Para Kᵖ el precondicionador sólo tiene sentido sobre el complemento del
    # espacio nulo (autovalores espurios ≈0 se excluyen).
    Dm = Diagonal(1.0 ./ sqrt.(diag(Me)))
    Mp = Dm * Me * Dm
    explica("Tabla III del paper: cond(Mᵖ) medido = $(round(cond(Mp), digits=2)) (esperado 19). ")
    @test isapprox(cond(Mp), 19.0, atol=0.5)

    Dk = Diagonal(1.0 ./ sqrt.(diag(Ke) .+ (diag(Ke) .== 0)))
    Kp = Dk * Ke * Dk
    λKp_nz = filter(x -> x > 1e-6, eigvals(Symmetric(Kp)))
    explica("cond(Kᵖ) sobre el complemento del núcleo espurio = $(round(maximum(λKp_nz) / minimum(λKp_nz), digits=2)) (esperado 35, Tabla III).")
    @test isapprox(maximum(λKp_nz) / minimum(λKp_nz), 35.0, atol=1.0)
end
