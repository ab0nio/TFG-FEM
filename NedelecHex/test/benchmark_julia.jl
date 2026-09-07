#!/usr/bin/env julia
# ===========================================================================
# Benchmark de rendimiento ensamblaje y solver en Julia
# ===========================================================================
#
# Etapas medidas : (coefficients, dof_numbering, assembly, solver)
# El assembly incluye la imposición de Dirichlet (dirichlet=true), igual que el
# assembly de MATLAB (getGlobalMatrices llama a imposeLocalDirichletMatrix).
# Metodología:
#   - BLAS (Basic Linear Algebra Subprograms) a 1 hilo.
#   - Se reportan por separado:
#       * TTFX (time-to-first-execution) que incluye la compilación JIT de Julia.
#       * mediana de Benchmark después.
#   - Memoria: bytes asignados en la etapa (allocations, medido por BenchmarkTools)
#     y memoria ocupada por las matrices dispersas resultantes (Base.summarysize).
#
# Mallas: las 5 exportadas de MATLAB (output/matlab_{nodes,elems}_N.csv) con matlab/exportMesh.m .
#
# Salida: output/results_julia.csv .
#
# Uso:  julia --project=. test/benchmark_julia.jl
# ===========================================================================

using NedelecHex
using BenchmarkTools, LinearAlgebra, SparseArrays, Printf

using NedelecHex.Mesh: load_matlab_mesh
using NedelecHex.GlobalDOFs: build_dof_map, dirichlet_dofs
using NedelecHex.Assembly: assemble_global_matrices, free_dofs
using NedelecHex.CoefficientGenerator: get_coefficients
using NedelecHex.Monomials: build_h2_monomials, build_h2_curl_monomials
using NedelecHex.HalfFilledCavity: analytic_eigenvalues, _lowest_physical_eigenvalues

# Un solo hilo BLAS.
BLAS.set_num_threads(1)

const N_MESHES = 5
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = normpath(joinpath(SCRIPT_DIR, "..", "output"))
const CSV_PATH = joinpath(OUTPUT_DIR, "results_julia.csv")

# Escribe una fila del CSV con el esquema común. Los campos no aplicables (p. ej.
# la memoria ocupada de una etapa que no produce matriz) se pasan como `nothing`
# y se escriben "NA".
#
# nnz_free / mem_free_bytes se refieren al bloque tras Dirichlet.

#funcion de formateo del CSV: si el valor es `nothing` → "NA", si no → string(valor)
function emit(io, malla, n_elem, n_dof, etapa,
    t_ttfx_ms, t_warm_ms, bytes_alloc, mem_occ_bytes, nnz_val,
    nnz_free=nothing, mem_free_bytes=nothing)
    na(x) = x === nothing ? "NA" : string(x) #los campos no aplicables se escriban como "NA" en el CSV
    @printf(io, "julia,%d,%d,%d,%s,%.6f,%.6f,%s,%s,%s,%s,%s\n",
        malla, n_elem, n_dof, etapa,
        t_ttfx_ms, t_warm_ms,
        na(bytes_alloc), na(mem_occ_bytes), na(nnz_val),
        na(nnz_free), na(mem_free_bytes))
end

function run_stages(io)
    println(io, "lenguaje,malla,n_elem,n_dof,etapa,t_ttfx_ms,t_warm_ms," *
                "bytes_alloc,mem_occupied_bytes,nnz,nnz_free,mem_free_bytes")

    # -------------------------------------------------------------------
    # Etapa 1: coeficientes (independiente de la malla → malla = 0)
    # -------------------------------------------------------------------
    println(stderr, "[julia] etapa coefficients …")
    t_ttfx = @elapsed get_coefficients()            # primera llamada: incluye JIT
    b = @benchmark get_coefficients() samples=200
    W, _ = get_coefficients()
    emit(io, 0, 0, 0, "coefficients",
        t_ttfx * 1e3, median(b).time / 1e6, median(b).memory, nothing, nothing)

    # Recursos compartidos por el ensamblaje (se construyen una vez; su coste es
    # parte de "coefficients", no del bucle de mallas).
    mono = build_h2_monomials()
    curl_mono = build_h2_curl_monomials()

    # Autovalores analíticos (banda física): sigma del solver eigs se centra aquí.
    k2_analytic = analytic_eigenvalues() .^ 2

    # -------------------------------------------------------------------
    # Etapas 2–4 por malla: dof_numbering, assembly (con Dirichlet), solver
    # -------------------------------------------------------------------
    for i in 1:N_MESHES
        # mallas de matlab, con sus etiquetas de material y de contorno
        mesh = load_matlab_mesh(OUTPUT_DIR, i)
        n_elem = size(mesh.elements, 1)
        println(stderr, "[julia] malla $i ($n_elem elem) …")

        # --- dof_numbering -------------------------------------------------
        t_ttfx = @elapsed build_dof_map(mesh)
        b = @benchmark build_dof_map($mesh) samples=100
        dof_map = build_dof_map(mesh)
        n_dof = dof_map.n_dofs
        emit(io, i, n_elem, n_dof, "dof_numbering",
            t_ttfx * 1e3, median(b).time / 1e6, median(b).memory, nothing, nothing)

        # --- assembly (incluye Dirichlet, como MATLAB) ---------------------
        # `dirichlet=true` impone las condiciones esenciales (PEC) POR ELEMENTO
        # dentro del ensamblaje, igual que MATLAB (+assembly/getGlobalMatrices.m
        # llama a imposeLocalDirichletMatrix). Así el assembly de Julia es
        # comparable 1:1 con el de MATLAB y la matriz "completa" devuelta ya
        # llega con los acoplamientos de contorno a cero.
        t_ttfx = @elapsed assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=true)
        b = @benchmark assemble_global_matrices($mesh, $dof_map, $W, $mono, $curl_mono; dirichlet=true) samples=50
        K, M = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=true)
        mem_occ = Base.summarysize(K) + Base.summarysize(M)

        # Kfree y Mfree para comparar con MATLAB
        free_stats = free_dofs(dof_map.n_dofs, dirichlet_dofs(mesh, dof_map))
        Kf_stats, Mf_stats = K[free_stats, free_stats], M[free_stats, free_stats]
        mem_free = Base.summarysize(Kf_stats) + Base.summarysize(Mf_stats)
        emit(io, i, n_elem, n_dof, "assembly",
            t_ttfx * 1e3, median(b).time / 1e6, median(b).memory, mem_occ, nnz(K),
            nnz(Kf_stats), mem_free)

        # --- solver (autovalores) ------------------------------------------
        # εᵣ de las etiquetas de la malla de MATLAB, para que sea el mismo problema
        Ks, Ms = assemble_global_matrices(mesh, dof_map, W, mono, curl_mono; dirichlet=true)
        free = free_dofs(dof_map.n_dofs, dirichlet_dofs(mesh, dof_map))
        Kf, Mf = Ks[free, free], Ms[free, free]
        # Dos muestras con @elapsed porque con benchmark me tardaba mucho
        solve() = _lowest_physical_eigenvalues(Kf, Mf, 6, k2_analytic)
        t_ttfx = @elapsed solve()      # incluye JIT
        t_warm = @elapsed solve()      # ya compilado
        emit(io, i, n_elem, n_dof, "solver",
            t_ttfx * 1e3, t_warm * 1e3, nothing, nothing, nothing)
    end
end

function main()
    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

    # Comprobar que las mallas exportadas existen ANTES de abrir el CSV, para no
    # dejar un fichero de resultados vacío si falta una malla.
    for i in 1:N_MESHES
        nf = joinpath(OUTPUT_DIR, "matlab_nodes_$i.csv")
        ef = joinpath(OUTPUT_DIR, "matlab_elems_$i.csv")
        if !isfile(nf) || !isfile(ef)
            error("Falta la malla $i ($nf). Genera los CSV con matlab/exportMesh.m.")
        end
    end

    # `open(...) do io` cierra el fichero aunque una etapa lance un error
    open(CSV_PATH, "w") do io
        run_stages(io)
    end
    println(stderr, "[julia] listo → $CSV_PATH")
end

main()
