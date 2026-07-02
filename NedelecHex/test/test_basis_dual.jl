@testset "BasisEval.jl — Propiedad dual de la base y curl" begin
    W, _      = NedelecHex.CoefficientGenerator.get_coefficients()
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    # ── Propiedad dual: gᵢ(Nⱼ) = δᵢⱼ ────────────────────────────────────────
    # A diferencia de test_G_matrix.jl (que sólo comprueba G·W=I algebraicamente),
    # aquí se evalúan las funciones de base reales vía eval_bases_ref y se les
    # aplican los funcionales DOF de DOFs.jl directamente, ejercitando el
    # camino completo BasisEval.jl → DOFs.jl de punta a punta.
    D = zeros(54, 54)
    for j in 1:54
        Nj  = pt -> NedelecHex.BasisEval.eval_bases_ref(W, mono, pt)[j, :]
        row = 1
        for edge_idx in 1:12
            D[row:row+1, j] = NedelecHex.DOFs.evaluate_edge_dof(Nj, edge_idx)
            row += 2
        end
        for face_idx in 1:6
            D[row:row+3, j] = NedelecHex.DOFs.evaluate_face_dof(Nj, face_idx)
            row += 4
        end
        D[row:row+5, j] = NedelecHex.DOFs.evaluate_volume_dof(Nj)
    end
    @test isapprox(D, I(54), atol=1e-8)

    # ── Curl de las funciones de base asociadas al volumen (49–54) ─────────
    # En el centro del elemento [0.5,0.5,0.5] el curl se anula por simetría,
    # así que se evalúa en un punto asimétrico.
    pt_asym = [0.2, 0.7, 0.9]
    curlN   = NedelecHex.BasisEval.eval_curl_bases_ref(W, curl_mono, pt_asym)
    for j in 49:54
        @test norm(curlN[j, :]) > 1e-6
    end

    # ── Curl de los modos irrotacionales (espurios) ──────────────────────────
    # El espacio nulo de Kᵉ = ∫ curlNᵢ·curlNⱼ corresponde a combinaciones
    # "gradiente" de las funciones de base cuyo curl es idénticamente nulo.
    # Se extrae un vector propio del autovalor más pequeño de Kᵉ y se
    # reconstruye el campo curl en varios puntos: debe ser ≈0 en todos ellos.
    nodes = NedelecHex.ReferenceHex.NODES
    Ke    = NedelecHex.Assembly.stiff_matrix_element(W, curl_mono, nodes)
    eg    = eigen(Symmetric(Ke))
    @test eg.values[1] < 1e-8         # autovalor nulo (espurio)
    x_null = eg.vectors[:, 1]

    for pt in ([0.5, 0.5, 0.5], [0.2, 0.7, 0.9], [0.6, 0.1, 0.4])
        curl_field = NedelecHex.BasisEval.eval_curl_bases_ref(W, curl_mono, pt)' * x_null
        @test norm(curl_field) < 1e-8
    end
end
