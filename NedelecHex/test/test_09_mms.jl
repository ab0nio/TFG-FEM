@testset "MMS — Exactitud en H₂ (Método de las Soluciones Manufacturadas)" begin
    # V_MMS = Σₗ mₗ es la suma de los 54 monomios de H₂, evaluada en
    # coordenadas físicas globales sobre el dominio [0,1]³. Restringido a
    # cualquier subelemento alineado con los ejes, sigue siendo un elemento
    # exacto del espacio local H₂ (la expansión binomial de un monomio afín
    # no aumenta su grado en ninguna variable), así que V_MMS es representable
    # de forma exacta en todo el espacio FEM global, sin importar la malla.
    #
    # Se manufactura el término de forzado fuerte O = curl(curl V_MMS) − k₀²V_MMS
    # (μᵣ=εᵣ=1) y se comprueba que, imponiendo como condición de contorno la
    # traza exacta de V_MMS en los DOFs de frontera y resolviendo el sistema
    # FEM con ese forzado, se recupera exactamente V_MMS: el sistema FEM es
    # consistente y el error debe ser ≈0 (nivel de precisión de máquina), no
    # sólo pequeño por refinamiento de malla.
    W, _          = NedelecHex.CoefficientGenerator.get_coefficients()
    mono          = NedelecHex.Monomials.build_h2_monomials()
    curl_mono     = NedelecHex.Monomials.build_h2_curl_monomials()
    curlcurl_mono = NedelecHex.Monomials.build_h2_curlcurl_monomials()

    V_MMS(pt) = sum(mono[l](pt) for l in 1:54)
    CCV(pt)   = sum(curlcurl_mono[l](pt) for l in 1:54)
    k0        = 1.0
    O(pt)     = CCV(pt) .- k0^2 .* V_MMS(pt)
    explica("MMS de exactitud: V_MMS = Σ de los 54 monomios pertenece a H₂ en cualquier malla alineada ⇒ el FEM debe recuperarla EXACTA (error ≈ precisión de máquina), no sólo 'pequeña'. El forzado O = ∇×∇×V − k₀²V se deriva analíticamente (curl-curl a mano).")

    # DOF exacto local de V_MMS para un elemento dado: se retrocede (pullback)
    # el campo físico a coordenadas de referencia del elemento — para un
    # hexaedro alineado con los ejes con Jacobiano diagonal J=diag(Δx,Δy,Δz),
    # el pullback covariante es V_ref = V_phys(x₀+J·ξ) · diag(J) (inverso del
    # Piola covariante N_phys = N_ref·J⁻¹ usado en piola_transform).
    function local_exact_dofs(nodes)
        J     = NedelecHex.BasisEval.compute_jacobian(nodes, [0.5, 0.5, 0.5])
        scale = [J[1,1], J[2,2], J[3,3]]
        x0    = nodes[1]
        V_ref = pt -> V_MMS(x0 .+ scale .* pt) .* scale

        d = zeros(54)
        row = 1
        for edge_idx in 1:12
            d[row:row+1] = NedelecHex.DOFs.evaluate_edge_dof(V_ref, edge_idx)
            row += 2
        end
        for face_idx in 1:6
            d[row:row+3] = NedelecHex.DOFs.evaluate_face_dof(V_ref, face_idx)
            row += 4
        end
        d[row:row+5] = NedelecHex.DOFs.evaluate_volume_dof(V_ref)
        return d
    end

    explica("Los DOFs exactos se calculan independientemente desde cada elemento que los toca y deben coincidir entre vecinos (V_MMS es H(curl)-conforme global); con la traza exacta impuesta en el contorno, el error relativo del sistema resuelto debe ser < 1e-9.")
    @printf("\n%-8s  %-10s  %s\n", "malla", "n_dofs", "‖x_fem − x_exacto‖ / ‖x_exacto‖")
    for n in (1, 2, 3)
        mesh    = NedelecHex.Mesh.box_mesh(nx=n, ny=n, nz=n)
        dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
        n_elem  = size(mesh.elements, 1)

        # DOFs exactos globales de V_MMS, calculados de forma independiente
        # desde cada elemento que toca cada DOF. Deben coincidir exactamente
        # entre elementos vecinos: es una comprobación adicional de que
        # V_MMS es realmente un campo H(curl) globalmente suave.
        x_exact = fill(NaN, dof_map.n_dofs)
        for e in 1:n_elem
            nidx   = mesh.elements[e, :]
            nodes  = [mesh.nodes[nidx[k], :] for k in 1:8]
            dexact = local_exact_dofs(nodes)
            dofs   = dof_map.connectivity[e, :]
            for i in 1:54
                g = dofs[i]
                if isnan(x_exact[g])
                    x_exact[g] = dexact[i]
                else
                    @test isapprox(x_exact[g], dexact[i], atol=1e-8)
                end
            end
        end

        K, M = NedelecHex.Assembly.assemble_global_matrices(mesh, dof_map, W, mono, curl_mono)
        F    = NedelecHex.Assembly.assemble_load_vector(mesh, dof_map, W, mono, O)
        A    = K - k0^2 .* M

        dir  = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
        free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)

        # Condensación estática: se imponen los DOFs de frontera al valor
        # exacto (no homogéneo, a diferencia de la condición PEC del
        # benchmark de la cavidad) y se resuelve sólo para los DOFs libres.
        rhs    = F[free] .- Matrix(A[free, dir]) * x_exact[dir]
        x_free = Matrix(A[free, free]) \ rhs

        x_fem        = copy(x_exact)
        x_fem[free] .= x_free

        ς = norm(x_fem .- x_exact) / norm(x_exact)
        @printf("%-8d  %-10d  %.3e\n", n, dof_map.n_dofs, ς)
        @test ς < 1e-9
    end
end
@testset "MMS — Tasa de convergencia (onda plana oblicua)" begin
    # Segundo test MMS: V_MMS = V_pol·e^{-jk₀(k̂·r)}, una onda plana transversal
    # (k̂·V_pol = 0) que satisface exactamente curl(curl V_MMS) = k₀²V_MMS
    # (identidad vectorial k̂×(k̂×V_pol) = -V_pol para V_pol⊥k̂, |k̂|=1), así que
    # el forzado fuerte O = curl curl V_MMS − k₀²V_MMS es idénticamente cero:
    # no hace falta ensamblar vector de carga, sólo imponer la traza compleja
    # exacta en los DOFs de frontera y resolver el sistema interior.
    #
    # A diferencia de T17, V_MMS NO es polinómica: no pertenece exactamente a
    # H₂, así que aparece un error de discretización genuino que debe decrecer
    # con la malla a la tasa teórica (2–3, según el enunciado) en vez de ser
    # ≈0 en cualquier malla.
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    k̂    = [1.0, 1.0, 1.0] / sqrt(3)
    Vpol = [1.0, -1.0, 0.0] / sqrt(2)
    @test isapprox(dot(k̂, Vpol), 0.0, atol=1e-14)
    k0 = 2.0
    explica("MMS de convergencia: onda plana transversal (k̂·V_pol = 0 ⇒ ∇×∇×V = k₀²V ⇒ forzado nulo). V no es polinómica, así que hay error de discretización genuino que debe caer con pendiente ≈ 2 en log-log, la tasa del elemento de 2º orden.")

    V_MMS(pt)      = Vpol .* exp(-1im * k0 * dot(k̂, pt))
    curl_V_MMS(pt) = (-1im * k0 * exp(-1im * k0 * dot(k̂, pt))) .* cross(k̂, Vpol)

    # Los funcionales DOF de arista/cara/volumen (DOFs.jl) sólo manejan campos
    # reales. Aquí se replican para campos complejos, con una salvedad crucial:
    # LinearAlgebra.dot conjuga su primer argumento para vectores complejos
    # (producto sesquilineal), lo que introduciría una conjugación espuria del
    # campo en cada DOF. Se usa sum(campo .* vector_real) en su lugar, que es
    # el producto bilineal (no conjugado) que requiere la definición del DOF.
    function complex_edge_dof(field, edge_idx)
        n1, n2 = NedelecHex.ReferenceHex.EDGES[edge_idx]
        p1 = Float64.(NedelecHex.ReferenceHex.NODES[n1])
        p2 = Float64.(NedelecHex.ReferenceHex.NODES[n2])
        Jedge = p2 - p1
        tau = Jedge / norm(Jedge)
        len = norm(Jedge)
        weights, points = NedelecHex.GaussQuadrature.gauss_segment(4)
        dofs = zeros(ComplexF64, 2)
        for qi in 1:2, k in eachindex(weights)
            t = points[k]
            pt3D = p1 + Jedge * t
            dofs[qi] += weights[k] * sum(field(pt3D) .* tau) * NedelecHex.DOFs.Q_EDGE[qi](t) * len
        end
        return dofs
    end
    function complex_face_dof(field, face_idx)
        weights, points = NedelecHex.GaussQuadrature.gauss_quadrilateral(4)
        dofs = zeros(ComplexF64, 4)
        for qi in 1:4, k in eachindex(weights)
            pt3D = NedelecHex.CoordinateTransforms.transform_2d_to_3d(points[k,:], face_idx)
            dofs[qi] += weights[k] * sum(field(pt3D) .* NedelecHex.DOFs.Q_FACE[face_idx][qi](pt3D))
        end
        return dofs
    end
    function complex_volume_dof(field)
        weights, points = NedelecHex.GaussQuadrature.gauss_hex(4)
        dofs = zeros(ComplexF64, 6)
        for qi in 1:6, k in eachindex(weights)
            pt = points[k,:]
            dofs[qi] += weights[k] * sum(field(pt) .* NedelecHex.DOFs.Q_VOL[qi](pt))
        end
        return dofs
    end

    # Pullback covariante V_ref = V_phys(x₀+J·ξ)·diag(J), igual que en T17.
    function local_exact_dofs(nodes)
        J     = NedelecHex.BasisEval.compute_jacobian(nodes, [0.5, 0.5, 0.5])
        scale = [J[1,1], J[2,2], J[3,3]]
        x0    = nodes[1]
        V_ref = pt -> V_MMS(x0 .+ scale .* pt) .* scale

        d = zeros(ComplexF64, 54)
        row = 1
        for edge_idx in 1:12
            d[row:row+1] = complex_edge_dof(V_ref, edge_idx)
            row += 2
        end
        for face_idx in 1:6
            d[row:row+3] = complex_face_dof(V_ref, face_idx)
            row += 4
        end
        d[row:row+5] = complex_volume_dof(V_ref)
        return d
    end

    ns     = [1, 2, 4, 8]
    hs     = 1.0 ./ ns
    ς_vals = zeros(length(ns))
    ς_curl_vals = zeros(length(ns))

    for (idx, n) in enumerate(ns)
        mesh    = NedelecHex.Mesh.box_mesh(nx=n, ny=n, nz=n)
        dof_map = NedelecHex.GlobalDOFs.build_dof_map(mesh)
        n_elem  = size(mesh.elements, 1)

        x_exact = fill(ComplexF64(NaN), dof_map.n_dofs)
        for e in 1:n_elem
            nidx   = mesh.elements[e, :]
            nodes  = [mesh.nodes[nidx[k], :] for k in 1:8]
            dexact = local_exact_dofs(nodes)
            dofs   = dof_map.connectivity[e, :]
            for i in 1:54
                g = dofs[i]
                isnan(x_exact[g]) && (x_exact[g] = dexact[i])
            end
        end

        K, M = NedelecHex.Assembly.assemble_global_matrices(mesh, dof_map, W, mono, curl_mono)
        A    = K - k0^2 .* M

        dir  = NedelecHex.GlobalDOFs.dirichlet_dofs(mesh, dof_map)
        free = NedelecHex.Assembly.free_dofs(dof_map.n_dofs, dir)

        rhs    = -(Matrix(A[free, dir]) * x_exact[dir])   # O=0 ⇒ sin vector de carga
        x_free = Matrix(A[free, free]) \ rhs
        x_fem  = copy(x_exact)
        x_fem[free] .= x_free

        # Error L² relativo del campo y de su curl, vía cuadratura de Gauss
        # sobre todos los elementos, comparando la reconstrucción FEM (Piola)
        # contra el campo/curl analítico exacto en cada punto físico.
        gaussW, gaussPts = NedelecHex.GaussQuadrature.gauss_hex(4)
        num_f, den_f, num_c, den_c = 0.0, 0.0, 0.0, 0.0
        for e in 1:n_elem
            nidx  = mesh.elements[e, :]
            nodes = [mesh.nodes[nidx[k], :] for k in 1:8]
            dofs  = dof_map.connectivity[e, :]
            xe    = x_fem[dofs]
            J     = NedelecHex.BasisEval.compute_jacobian(nodes, [0.5, 0.5, 0.5])
            detJ  = abs(det(J))
            x0    = nodes[1]
            scale = [J[1,1], J[2,2], J[3,3]]
            for g in eachindex(gaussW)
                pt     = gaussPts[g, :]
                x_phys = x0 .+ scale .* pt

                Nref  = NedelecHex.BasisEval.eval_bases_ref(W, mono, pt)
                Nphys = NedelecHex.BasisEval.piola_transform(Nref, J)
                E_fem = Nphys' * xe
                E_ex  = V_MMS(x_phys)
                num_f += gaussW[g] * detJ * sum(abs2, E_fem .- E_ex)
                den_f += gaussW[g] * detJ * sum(abs2, E_ex)

                cNref  = NedelecHex.BasisEval.eval_curl_bases_ref(W, curl_mono, pt)
                cNphys = NedelecHex.BasisEval.curl_piola_transform(cNref, J)
                curlE_fem = cNphys' * xe
                curlE_ex  = curl_V_MMS(x_phys)
                num_c += gaussW[g] * detJ * sum(abs2, curlE_fem .- curlE_ex)
                den_c += gaussW[g] * detJ * sum(abs2, curlE_ex)
            end
        end
        ς_vals[idx]      = sqrt(num_f / den_f)
        ς_curl_vals[idx] = sqrt(num_c / den_c)
    end

    @printf("\n%-8s  %-12s  %-12s  %-12s\n", "n", "h", "ς (campo)", "ς_curl")
    for i in eachindex(ns)
        @printf("%-8d  %-12.4f  %-12.3e  %-12.3e\n", ns[i], hs[i], ς_vals[i], ς_curl_vals[i])
    end

    slope(y) = begin
        x  = log.(hs)
        yy = log.(y)
        x̄, ȳ = sum(x)/length(x), sum(yy)/length(yy)
        sum((x .- x̄) .* (yy .- ȳ)) / sum((x .- x̄).^2)
    end
    slope_field = slope(ς_vals)
    slope_curl  = slope(ς_curl_vals)
    @printf("\nPendiente log-log: ς=%.2f  ς_curl=%.2f  (esperada entre 2 y 3)\n", slope_field, slope_curl)

    @test issorted(ς_vals, rev=true)
    @test issorted(ς_curl_vals, rev=true)
    @test 1.8 <= slope_field <= 3.2
    @test 1.8 <= slope_curl  <= 3.2
end
