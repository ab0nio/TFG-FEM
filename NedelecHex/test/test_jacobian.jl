@testset "BasisEval.jl — Jacobiano del elemento de referencia" begin
    nodes = NedelecHex.ReferenceHex.NODES

    # Para el hexaedro de referencia [0,1]³ con nodos exactamente en NODES,
    # el mapeo isoparamétrico es la identidad, así que J debe ser I₃ en
    # cualquier punto de referencia.
    for pt in ([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.25, 0.75, 0.5])
        J = NedelecHex.BasisEval.compute_jacobian(nodes, pt)
        @test norm(J - I) < 1e-12
        @test isapprox(det(J), 1.0, atol=1e-12)
        @test norm(inv(J)' - I) < 1e-12
    end

    # Para un hexaedro escalado [0,a]×[0,b]×[0,c], J debe ser diag(a,b,c),
    # det(J) = a·b·c y J⁻ᵀ = diag(1/a,1/b,1/c) (J diagonal ⇒ J⁻ᵀ = J⁻¹).
    a, b, c = 2.0, 3.0, 0.5
    scaled_nodes = [n .* [a, b, c] for n in NedelecHex.ReferenceHex.NODES]
    for pt in ([0.5, 0.5, 0.5], [0.1, 0.9, 0.3])
        J_scaled = NedelecHex.BasisEval.compute_jacobian(scaled_nodes, pt)
        @test norm(J_scaled - Diagonal([a, b, c])) < 1e-12
        @test isapprox(det(J_scaled), a * b * c, atol=1e-10)
        @test norm(inv(J_scaled)' - Diagonal([1/a, 1/b, 1/c])) < 1e-10
    end
end
