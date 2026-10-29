# Genera las etiquetas de material de las mallas del paper a partir de la
# geometría: la losa es el semiespacio x > (x0+x1)/2 (εᵣ=2), el resto vacío (εᵣ=1).
# Escribe, por malla i, matlab_solids_i.csv, matlab_elemsolid_i.csv y
# matlab_eps_i.csv en test/data/meshes/.
# Uso: julia --project=. scripts/generate_material_tags.jl

using DelimitedFiles

const MESH_DIR = joinpath(@__DIR__, "..", "test", "data", "meshes")

function generate(i::Int)
    nodes = readdlm(joinpath(MESH_DIR, "matlab_nodes_$i.csv"), ',', Float64)
    elems = readdlm(joinpath(MESH_DIR, "matlab_elems_$i.csv"), ',', Int)
    n_elem = size(elems, 1)

    x0, x1 = extrema(@view nodes[:, 1])
    mid = (x0 + x1) / 2

    solid = Vector{Int}(undef, n_elem)
    for e in 1:n_elem
        x_c = sum(@view nodes[elems[e, :], 1]) / 8
        solid[e] = x_c > mid ? 2 : 1
    end

    open(joinpath(MESH_DIR, "matlab_solids_$i.csv"), "w") do io
        println(io, "1,1.0,1.0,vacio")
        println(io, "2,2.0,1.0,dielectrico")
    end
    open(joinpath(MESH_DIR, "matlab_elemsolid_$i.csv"), "w") do io
        for s in solid
            println(io, s)
        end
    end
    open(joinpath(MESH_DIR, "matlab_eps_$i.csv"), "w") do io
        for s in solid
            println(io, s == 2 ? "2.0" : "1.0")
        end
    end

    n2 = count(==(2), solid)
    println("malla $i: n_elem=$n_elem, εᵣ=2 → $n2, mitad=$(n_elem ÷ 2), " *
        (n2 == n_elem ÷ 2 ? "OK" : "MISMATCH"))
    return n2 == n_elem ÷ 2
end

ok = all(generate(i) for i in 1:5)
ok || error("el recuento de elementos con εᵣ=2 no es la mitad en alguna malla")
