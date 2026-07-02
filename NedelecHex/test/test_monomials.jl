@testset "Monomials.jl — Los 54 monomios vectoriales de H₂" begin
    mono      = NedelecHex.Monomials.build_h2_monomials()
    curl_mono = NedelecHex.Monomials.build_h2_curl_monomials()

    @test length(mono) == 54
    @test length(curl_mono) == 54

    # ── 18 monomios por componente (ξ̂: 1–18, η̂: 19–36, ζ̂: 37–54) ────────────
    samples = ([0.2, 0.6, 0.9], [0.7, 0.1, 0.4], [0.35, 0.85, 0.15])
    for pt in samples
        for j in 1:18
            v = mono[j](pt)
            @test v[2] == 0.0 && v[3] == 0.0
        end
        for j in 19:36
            v = mono[j](pt)
            @test v[1] == 0.0 && v[3] == 0.0
        end
        for j in 37:54
            v = mono[j](pt)
            @test v[1] == 0.0 && v[2] == 0.0
        end
    end

    # ── Restricción de Nédélec sobre el grado polinómico ────────────────────
    # n-ésima diferencia finita de f respecto a la variable x en x0, paso h:
    # es exactamente 0 (dentro de la precisión de máquina) si deg(f) < n.
    findiff(f, x0, h, n) =
        sum((-1)^(n - k) * binomial(n, k) * f(x0 + k * h) for k in 0:n)

    # component: función escalar de (ξ,η,ζ), axis: 1=ξ,2=η,3=ζ
    function scalar_along(component, axis, fixed)
        return t -> begin
            pt = copy(fixed)
            pt[axis] = t
            component(pt)
        end
    end

    fixed = [0.37, 0.61, 0.83]
    h = 0.05

    # Bloque ξ̂ (1–18): grado ≤1 en ξ, ≤2 en η y ζ
    for j in 1:18
        g = pt -> mono[j](pt)[1]
        @test isapprox(findiff(scalar_along(g, 1, fixed), fixed[1], h, 2), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 2, fixed), fixed[2], h, 3), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 3, fixed), fixed[3], h, 3), 0.0, atol=1e-10)
    end

    # Bloque η̂ (19–36): grado ≤2 en ξ, ≤1 en η, ≤2 en ζ
    for j in 19:36
        g = pt -> mono[j](pt)[2]
        @test isapprox(findiff(scalar_along(g, 1, fixed), fixed[1], h, 3), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 2, fixed), fixed[2], h, 2), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 3, fixed), fixed[3], h, 3), 0.0, atol=1e-10)
    end

    # Bloque ζ̂ (37–54): grado ≤2 en ξ, ≤2 en η, ≤1 en ζ
    for j in 37:54
        g = pt -> mono[j](pt)[3]
        @test isapprox(findiff(scalar_along(g, 1, fixed), fixed[1], h, 3), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 2, fixed), fixed[2], h, 3), 0.0, atol=1e-10)
        @test isapprox(findiff(scalar_along(g, 3, fixed), fixed[3], h, 2), 0.0, atol=1e-10)
    end

    # ── Consistencia curl analítico vs. diferencias finitas centradas ─────────
    function fd_curl(f, pt; hh=1e-4)
        function pd(idx, axis)
            p1 = copy(pt); p1[axis] += hh
            p2 = copy(pt); p2[axis] -= hh
            return (f(p1)[idx] - f(p2)[idx]) / (2hh)
        end
        return [pd(3, 2) - pd(2, 3), pd(1, 3) - pd(3, 1), pd(2, 1) - pd(1, 2)]
    end

    pt_test = [0.3, 0.6, 0.8]
    for j in 1:54
        c_analytic = curl_mono[j](pt_test)
        c_fd       = fd_curl(mono[j], pt_test)
        @test isapprox(c_analytic, c_fd, atol=1e-6)
    end
end
