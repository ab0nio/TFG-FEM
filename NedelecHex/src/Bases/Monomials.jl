module Monomials

export build_h2_monomials

function build_h2_monomials()
    m = Vector{Function}()

    # Bloque X: componente [*, 0, 0]
    # MATLAB getMixedPolySpace([x,y,z], order1=1, order2=2, order3=2)
    # Bucle: ii=kz+1 (exterior), jj=ky+1 (medio), kk=kx+1 (interior)
    # Condición: (ii+jj+kk-3) <= (1+2+2) → siempre se cumple para estos rangos
    for kz in 0:2, ky in 0:2, kx in 0:1
        let kx = kx, ky = ky, kz = kz
            push!(m, p -> [p[1]^kx * p[2]^ky * p[3]^kz, 0.0, 0.0])
        end
    end

    # Bloque Y: componente [0, *, 0]
    # MATLAB getMixedPolySpace([x,y,z], order1=2, order2=1, order3=2)
    # Bucle: ii=kz+1 (exterior), jj=ky+1 (medio), kk=kx+1 (interior)
    for kz in 0:2, ky in 0:1, kx in 0:2
        let kx = kx, ky = ky, kz = kz
            push!(m, p -> [0.0, p[1]^kx * p[2]^ky * p[3]^kz, 0.0])
        end
    end

    # Bloque Z: componente [0, 0, *]
    # MATLAB getMixedPolySpace([x,y,z], order1=2, order2=2, order3=1)
    # Bucle: ii=kz+1 (exterior), jj=ky+1 (medio), kk=kx+1 (interior)
    for kz in 0:1, ky in 0:2, kx in 0:2
        let kx = kx, ky = ky, kz = kz
            push!(m, p -> [0.0, 0.0, p[1]^kx * p[2]^ky * p[3]^kz])
        end
    end

    return m
end

end # module Monomials