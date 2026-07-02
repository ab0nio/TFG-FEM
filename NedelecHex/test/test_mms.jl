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
