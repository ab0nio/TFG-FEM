# Monomios de H₂ = Q₁,₂,₂ × Q₂,₁,₂ × Q₂,₂,₁
# Cada monomio se representa por la terna de exponentes (ξ, η, ζ)
# de su única componente no nula.

# Producto tensorial de un bloque Q_{gξ,gη,gζ}: todos los monomios
# ξ^a η^b ζ^c con 0 ≤ a ≤ gξ, 0 ≤ b ≤ gη, 0 ≤ c ≤ gζ
function bloque_Q(gξ, gη, gζ)
    return [(a, b, c) for a in 0:gξ, b in 0:gη, c in 0:gζ][:]
end

# Los tres bloques de H₂ según la restricción de orden mixto
H₂_ξ = bloque_Q(1, 2, 2)   # componente ξ : Q₁,₂,₂  → 18 monomios
H₂_η = bloque_Q(2, 1, 2)   # componente η : Q₂,₁,₂  → 18 monomios
H₂_ζ = bloque_Q(2, 2, 1)   # componente ζ : Q₂,₂,₁  → 18 monomios

# Base completa: 54 monomios vectoriales
H₂ = (ξ=H₂_ξ, η=H₂_η, ζ=H₂_ζ)

# Comprobación de la dimensión
@assert length(H₂.ξ) == 18 && length(H₂.η) == 18 && length(H₂.ζ) == 18
println("OK: 3 bloques × 18 monomios = ", length(H₂.ξ) + length(H₂.η) + length(H₂.ζ))
