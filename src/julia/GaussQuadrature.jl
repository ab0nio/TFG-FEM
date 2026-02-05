"""
Estructura para almacenar raíces y pesos de Gauss-Legendre
"""
struct GaussLegendreData
    roots::Vector{Float64}
    weights::Vector{Float64}
end

"""
tabla tabulada de raíces y pesos para n = 1 a 10
"""
const LEGENDRE_TABLE = Dict{Int,GaussLegendreData}(
    1 => GaussLegendreData(
        [0.0],
        [2.0]
    ),
    2 => GaussLegendreData(
        [-0.5773502691896257, 0.5773502691896257],
        [1.0, 1.0]
    ),
    3 => GaussLegendreData(
        [-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]
    ),
    4 => GaussLegendreData(
        [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]
    ),
    5 => GaussLegendreData(
        [-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891]
    ),
    6 => GaussLegendreData(
        [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521],
        [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704]
    ),
    7 => GaussLegendreData(
        [-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585],
        [0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697]
    ),
    8 => GaussLegendreData(
        [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363],
        [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]
    ),
    9 => GaussLegendreData(
        [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261],
        [0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, 0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.0812743883615744]
    ),
    10 => GaussLegendreData(
        [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717],
        [0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881]
    )
)

function gauss_legendre_integrate(f::Function, n::Int; a::Real=-1.0, b::Real=1.0)
    # Validar entrada
    if n < 1 || n > 10
        error("n debe estar entre 1 y 10.")
    end

    # Obtener raíces y pesos tabulados
    gauss_data = LEGENDRE_TABLE[n]
    xi_standard = gauss_data.roots
    wi_standard = gauss_data.weights

    # Transformación de intervalo 
    if (a, b) != (-1.0, 1.0)
        scale_factor = (b - a) / 2.0
        shift = (b + a) / 2.0

        xi_transformed = scale_factor .* xi_standard .+ shift

        wi_scaled = scale_factor .* wi_standard
    else
        xi_transformed = xi_standard
        wi_scaled = wi_standard
    end


    integral = sum(wi_scaled .* f.(xi_transformed))

    return integral
end


res1 = gauss_legendre_integrate(x -> x^2, 1, a=5, b=10)

println(" ∫x² [5,10] Debe salir aproximadamente 291.6666666666667")
println(res1)

res2 = gauss_legendre_integrate(x -> sin(x)^2, 8, a=0, b=π)

println(" ∫sin²(x) [0,π] Debe salir aproximadamente 1.5707963267948966")
println(res2)

#lo que no se es como se sabe cuantos puntos hacen falta en funciones no polinomicas

res3 = gauss_legendre_integrate(x -> exp(x), 10, a=0, b=π)

println(" ∫e^x [0,π] Debe salir aproximadamente 22.14069263277927")
println(res3)
