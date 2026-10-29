# =============================================================================
#  sumsq.jl — Evaluación de una función de base hexaédrica de Nédélec
#  (componente ξ de H2) sobre una nube de puntos del hexaedro de referencia.
#
#  Julia permite usar directamente caracteres griegos (ξ, η, ζ) como
#  identificadores, igual que en la notación matemática del TFG.
#
#  Ejecutar:  julia sumsq.jl
# =============================================================================

# Función de base: monomio ξ·η²·ζ² (uno de los 54 de H2)
basis(ξ, η, ζ) = ξ * η^2 * ζ^2

function evaluate(n::Int)
    total = 0.0
    # Recorremos una nube n×n×n de puntos del hexaedro de referencia [0,1]^3
    for i in 1:n
        ξ = i / n
        for j in 1:n
            η = j / n
            for k in 1:n
                ζ = k / n
                total += basis(ξ, η, ζ)
            end
        end
    end
    return total
end

# Calentamiento (compilación JIT) para medir solo el cómputo
evaluate(10)

n = 200
t = @elapsed result = evaluate(n)
println("Julia   n=$n  resultado=$(result)  tiempo=$(round(t*1000, digits=2)) ms")
