module Monomials

export build_h2_monomials, build_h2_curl_monomials, build_h2_curlcurl_monomials

@inline _dpow(p, k) = k == 0 ? 0.0 : k * p^(k - 1)
@inline _d2pow(p, k) = k < 2 ? 0.0 : k * (k - 1) * p^(k - 2)

# 54 monomios vectoriales del espacio H₂ (ec. 4 del paper).
# Producto tensorial Q(1,2,2) × Q(2,1,2) × Q(2,2,1), en ese orden.
function build_h2_monomials()
    m = Vector{Function}()
    for kz in 0:2, ky in 0:2, kx in 0:1
        let kx=kx, ky=ky, kz=kz
            push!(m, p -> [p[1]^kx * p[2]^ky * p[3]^kz, 0.0, 0.0])
        end
    end
    for kz in 0:2, ky in 0:1, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(m, p -> [0.0, p[1]^kx * p[2]^ky * p[3]^kz, 0.0])
        end
    end
    for kz in 0:1, ky in 0:2, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(m, p -> [0.0, 0.0, p[1]^kx * p[2]^ky * p[3]^kz])
        end
    end
    return m
end

# Curl analítico de los 54 monomios. Mismo orden que build_h2_monomials.
# [f,0,0]→[0, ∂f/∂z, −∂f/∂y]   [0,g,0]→[−∂g/∂z, 0, ∂g/∂x]   [0,0,h]→[∂h/∂y, −∂h/∂x, 0]
function build_h2_curl_monomials()
    c = Vector{Function}()
    for kz in 0:2, ky in 0:2, kx in 0:1
        let kx=kx, ky=ky, kz=kz
            push!(c, p -> [0.0,
                            _dpow(p[3],kz)*p[1]^kx*p[2]^ky,
                           -_dpow(p[2],ky)*p[1]^kx*p[3]^kz])
        end
    end
    for kz in 0:2, ky in 0:1, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(c, p -> [-_dpow(p[3],kz)*p[1]^kx*p[2]^ky,
                            0.0,
                            _dpow(p[1],kx)*p[2]^ky*p[3]^kz])
        end
    end
    for kz in 0:1, ky in 0:2, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(c, p -> [ _dpow(p[2],ky)*p[1]^kx*p[3]^kz,
                            -_dpow(p[1],kx)*p[2]^ky*p[3]^kz,
                             0.0])
        end
    end
    return c
end


# Curl del curl analítico de los 54 monomios (necesario para el término de
# forzado O = curl(μᵣ⁻¹curl V) − k₀²εᵣV del método de las soluciones
# manufacturadas). Se deriva vía la identidad curl curl F = grad(div F) − ∇²F,
# aplicada a cada uno de los tres bloques de componente única. Mismo orden que
# build_h2_monomials.
function build_h2_curlcurl_monomials()
    cc = Vector{Function}()
    for kz in 0:2, ky in 0:2, kx in 0:1
        let kx=kx, ky=ky, kz=kz
            push!(cc, p -> [-_d2pow(p[2],ky)*p[1]^kx*p[3]^kz - _d2pow(p[3],kz)*p[1]^kx*p[2]^ky,
                              _dpow(p[1],kx)*_dpow(p[2],ky)*p[3]^kz,
                              _dpow(p[1],kx)*p[2]^ky*_dpow(p[3],kz)])
        end
    end
    for kz in 0:2, ky in 0:1, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(cc, p -> [ _dpow(p[1],kx)*_dpow(p[2],ky)*p[3]^kz,
                             -_d2pow(p[1],kx)*p[2]^ky*p[3]^kz - p[1]^kx*p[2]^ky*_d2pow(p[3],kz),
                              p[1]^kx*_dpow(p[2],ky)*_dpow(p[3],kz)])
        end
    end
    for kz in 0:1, ky in 0:2, kx in 0:2
        let kx=kx, ky=ky, kz=kz
            push!(cc, p -> [ _dpow(p[1],kx)*p[2]^ky*_dpow(p[3],kz),
                              p[1]^kx*_dpow(p[2],ky)*_dpow(p[3],kz),
                             -_d2pow(p[1],kx)*p[2]^ky*p[3]^kz - p[1]^kx*_d2pow(p[2],ky)*p[3]^kz])
        end
    end
    return cc
end

end # module Monomials
