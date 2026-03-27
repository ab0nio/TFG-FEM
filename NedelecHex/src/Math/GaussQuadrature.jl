module GaussQuadrature

export gauss_segment, gauss_quadrilateral, gauss_hex
#no usar valores tabulados
function gauss_segment(order::Int)
    if order <= 1
        weights = [1.0]
        points = [0.555555555555555]
    elseif order <= 3
        weights = [0.5, 0.5]
        points = [0.788675134594813, 0.211324865405187]
    elseif order <= 5
        weights = [0.277777777777778, 0.444444444444444, 0.277777777777778]
        points = [0.112701665379259, 0.500000000000000, 0.887298334620741]
    else
        error("Order not coded")
    end
    return weights, points
end

function gauss_quadrilateral(order::Int)
    w1d, p1d = gauss_segment(order)
    n = length(w1d)
    N = n^2

    weights = Vector{Float64}(undef, N)
    points = Matrix{Float64}(undef, N, 2)

    for i in 1:n
        idx = (i-1)*n+1:i*n
        weights[idx] = w1d[i] .* w1d
        points[idx, 1] = p1d
        points[idx, 2] .= p1d[i]
    end

    return weights, points
end

function gauss_hex(order::Int)
    w1d, p1d = gauss_segment(order)
    w2d, p2d = gauss_quadrilateral(order)

    n1d = length(w1d)
    n2d = length(w2d)
    N = n1d * n2d

    weights = Vector{Float64}(undef, N)
    points = Matrix{Float64}(undef, N, 3)

    for i in 1:n1d
        idx = (i-1)*n2d+1:i*n2d
        weights[idx] = w1d[i] .* w2d
        points[idx, 1:2] = p2d
        points[idx, 3] .= p1d[i]
    end

    return weights, points
end

end  # module GaussQuadrature