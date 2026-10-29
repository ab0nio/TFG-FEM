module BasisEval

using LinearAlgebra
using Printf

export eval_bases_ref, eval_curl_bases_ref,
    compute_jacobian, piola_transform, curl_piola_transform,
    save_bases_eval

# N = W' · M,  M[j,:] = m_j(pt).  Devuelve 54×3 (fila i = N_i(pt)).
function eval_bases_ref(W, monomials, pt)
    M = zeros(length(monomials), 3)
    for j in eachindex(monomials)
        M[j, :] = monomials[j](pt)
    end
    return W' * M
end

# Igual que eval_bases_ref pero con los curl de los monomios.
function eval_curl_bases_ref(W, curl_monomials, pt)
    C = zeros(length(curl_monomials), 3)
    for j in eachindex(curl_monomials)
        C[j, :] = curl_monomials[j](pt)
    end
    return W' * C
end

# J[i,j] = ∑_k (∂N_k/∂ξ_i)(pt) · x_k[j]
# Vale para cualquier hexaedro de 8 nodos, distorsionado o no; J es constante
# para hexaedros alineados con los ejes.
function compute_jacobian(nodes, pt)
    ξ, η, ζ = pt
    dNdξ = [-(1-η)*(1-ζ), (1-η)*(1-ζ), η*(1-ζ), -η*(1-ζ),
        -(1-η)*ζ, (1-η)*ζ, η*ζ, -η*ζ]
    dNdη = [-(1-ξ)*(1-ζ), -ξ*(1-ζ), ξ*(1-ζ), (1-ξ)*(1-ζ),
        -(1-ξ)*ζ, -ξ*ζ, ξ*ζ, (1-ξ)*ζ]
    dNdζ = [-(1-ξ)*(1-η), -ξ*(1-η), -ξ*η, -(1-ξ)*η,
        (1-ξ)*(1-η), ξ*(1-η), ξ*η, (1-ξ)*η]
    J = zeros(3, 3)
    for k in 1:8
        J[1, :] .+= dNdξ[k] * nodes[k]
        J[2, :] .+= dNdη[k] * nodes[k]
        J[3, :] .+= dNdζ[k] * nodes[k]
    end
    return J
end

# Piola covariante H(curl): jacobianMatrix\bases.'
function piola_transform(N_ref, J)
    return (J \ N_ref')'
end

# Piola del curl H(curl): curlBases * jacobianMatrix / det
function curl_piola_transform(curlN_ref, J)
    return curlN_ref * J / det(J)
end

function save_bases_eval(W, monomials, nodes, output_dir; pt=[0.5, 0.5, 0.5])
    N_ref = eval_bases_ref(W, monomials, pt)
    J = compute_jacobian(nodes, pt)
    N_phys = piola_transform(N_ref, J)
    mkpath(output_dir)
    open(joinpath(output_dir, "bases_eval_julia.txt"), "w") do f
        for i in 1:54
            @printf(f, "%.15f\t%.15f\t%.15f\n", N_phys[i, 1], N_phys[i, 2], N_phys[i, 3])
        end
    end
    println("Bases evaluadas guardadas en: $(joinpath(output_dir, "bases_eval_julia.txt"))")
end

end # module BasisEval
