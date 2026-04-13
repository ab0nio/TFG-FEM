module CoefficientGenerator

using LinearAlgebra
using BenchmarkTools: @benchmark, median, mean, minimum, maximum
using ..Monomials: build_h2_monomials
using ..DOFs: evaluate_edge_dof, evaluate_face_dof, evaluate_volume_dof
using Printf

export get_coefficients, run_all_tests

function get_coefficients()
    monomials = build_h2_monomials()
    N = length(monomials)
    G = zeros(N, N)

    for j in 1:N
        mono = monomials[j]
        row = 1

        for edge_idx in 1:12
            dofs = evaluate_edge_dof(mono, edge_idx)
            G[row:row+1, j] = dofs
            row += 2
        end

        for face_idx in 1:6
            dofs = evaluate_face_dof(mono, face_idx)
            G[row:row+3, j] = dofs
            row += 4
        end

        dofs = evaluate_volume_dof(mono)
        G[row:row+5, j] = dofs
        row += 6
    end

    W = G \ I(N)
    return W, G
end

function run_all_tests(; verbose::Bool=true)
    W, G = get_coefficients()

    # Validación mínima
    err = norm(W * G - I(54))
    @printf("W·G - I  error: %.2e  %s\n\n", err, err < 1e-10 ? "✓" : "✗")
    passed = err < 1e-10

    if verbose
        W_int = round.(Int, W)
        dof_labels = vcat(
            ["a$i" for i in 1:18],
            ["b$i" for i in 1:18],
            ["c$i" for i in 1:18]
        )
        sep = repeat("----|", 55)

        # Imprime por pantalla
        @printf("%-4s", "Ni")
        for j in 1:54
            @printf("|%4d", j)
        end
        println()
        println(sep)
        for i in 1:54
            @printf("%-4s", dof_labels[i])
            for j in 1:54
                @printf("|%4d", W_int[i, j])
            end
            println()
        end
        println(sep)

        # Guarda en fichero
        output_dir = joinpath(@__DIR__, "..", "..", "output")
        mkpath(output_dir)
        open(joinpath(output_dir, "tabla_W.txt"), "w") do f
            @printf(f, "%-4s", "Ni")
            for j in 1:54
                @printf(f, "|%4d", j)
            end
            println(f)
            println(f, sep)
            for i in 1:54
                @printf(f, "%-4s", dof_labels[i])
                for j in 1:54
                    @printf(f, "|%4d", W_int[i, j])
                end
                println(f)
            end
            println(f, sep)
        end
        println("Tabla W guardada en: $(joinpath(output_dir, "tabla_W.txt"))")
    end

    # Benchmark get_coefficients() and append results to the file
    println("\nInicio del benchmark:")
    bench = @benchmark get_coefficients() samples = 50 evals = 1
    t_min = minimum(bench).time / 1e6   # ns → ms
    t_median = median(bench).time / 1e6
    t_mean = mean(bench).time / 1e6
    t_max = maximum(bench).time / 1e6
    allocs = median(bench).allocs
    mem_kb = median(bench).memory / 1024

    @printf("Benchmark (n=%d muestras):\n", length(bench.times))
    @printf("  min:    %.3f ms\n", t_min)
    @printf("  median: %.3f ms\n", t_median)
    @printf("  mean:   %.3f ms\n", t_mean)
    @printf("  max:    %.3f ms\n", t_max)
    @printf("  allocs: %d  memory: %.1f KiB\n", allocs, mem_kb)

    output_dir = joinpath(@__DIR__, "..", "..", "output")
    open(joinpath(output_dir, "tabla_W.txt"), "a") do f
        println(f)
        println(f, "=== Benchmark: get_coefficients() (Julia) ===")
        @printf(f, "Samples: %d\n", length(bench.times))
        @printf(f, "  min:    %.3f ms\n", t_min)
        @printf(f, "  median: %.3f ms\n", t_median)
        @printf(f, "  mean:   %.3f ms\n", t_mean)
        @printf(f, "  max:    %.3f ms\n", t_max)
        @printf(f, "  allocs: %d  memory: %.1f KiB\n", allocs, mem_kb)
    end
    println("Benchmark guardado en: $(joinpath(output_dir, "tabla_W.txt"))")

    return passed, W, G
end

end  # module CoefficientGenerator
