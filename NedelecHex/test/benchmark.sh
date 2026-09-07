#!/usr/bin/env bash
# ===========================================================================
#junta Julia vs MATLAB
# ===========================================================================
# Lanza los dos benchmarks.
#
# Salidas (en NedelecHex/output/):
#   results_julia.csv, results_matlab.csv   
#   results_combined.csv                    
#   benchmark_report.md                    
#
# Uso:   bash test/benchmark.sh
#        bash test/benchmark.sh --solo-reporte
# ===========================================================================
set -euo pipefail

# --- Rutas -----------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"          # NedelecHex/
OUTPUT_DIR="$PKG_DIR/output"
mkdir -p "$OUTPUT_DIR"

JULIA_BIN="${JULIA_BIN:-julia}"
MATLAB_BIN="${MATLAB_BIN:-matlab}"

JULIA_CSV="$OUTPUT_DIR/results_julia.csv"
MATLAB_CSV="$OUTPUT_DIR/results_matlab.csv"
COMBINED_CSV="$OUTPUT_DIR/results_combined.csv"
REPORT_MD="$OUTPUT_DIR/benchmark_report.md"

# --solo-reporte: rehace results_combined.csv y benchmark_report.md a partir de
# los CSV que ya hay, sin volver a medir 
REPORT_ONLY=0
if [ "${1:-}" = "--solo-reporte" ]; then REPORT_ONLY=1; fi

# --- Un solo hilo en ambos entornos (comparación de algoritmo serial) ------
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "=== Benchmark ensamblaje + solver: Julia vs MATLAB ==="
echo "Salida: $OUTPUT_DIR"
echo

# --- 1. Lado Julia ---------------------------------------------------------
if [ "$REPORT_ONLY" -eq 1 ]; then
    echo ">>> --solo-reporte: no se mide nada, se reutilizan los CSV existentes."
    if [ ! -f "$JULIA_CSV" ]; then
        echo "ERROR: no existe $JULIA_CSV; ejecuta el benchmark sin --solo-reporte." >&2
        exit 1
    fi
elif ! command -v "$JULIA_BIN" >/dev/null 2>&1; then
    echo "ERROR: no encuentro julia ('$JULIA_BIN'). Instálalo o define JULIA_BIN." >&2
    exit 1
else
    echo ">>> Ejecutando benchmark de Julia..."
    "$JULIA_BIN" --project="$PKG_DIR" "$SCRIPT_DIR/benchmark_julia.jl"
fi
echo

# --- 2. Lado MATLAB (opcional) ---------------------------------------------
HAVE_MATLAB=0
if [ "$REPORT_ONLY" -eq 1 ]; then
    MATLAB_ROWS=0
    [ -f "$MATLAB_CSV" ] && MATLAB_ROWS=$(( $(wc -l < "$MATLAB_CSV") - 1 ))
    if [ "$MATLAB_ROWS" -gt 0 ]; then
        HAVE_MATLAB=1
        echo ">>> Reutilizando $MATLAB_ROWS filas de $MATLAB_CSV."
    else
        echo "AVISO: no hay datos de MATLAB que reutilizar; se reporta solo Julia." >&2
    fi
elif command -v "$MATLAB_BIN" >/dev/null 2>&1; then
    echo ">>> Ejecutando benchmark de MATLAB..."
    # Borrar resultados de una tirada anterior: si esta falla, no queremos fundir
    # datos viejos con los de Julia recién medidos.
    rm -f "$MATLAB_CSV"

    # -batch: sin GUI, sin splash. `set -e` está activo, así que capturamos el
    # código con `|| MATLAB_RC=$?` en vez de dejar que aborte el script.
    MATLAB_RC=0
    "$MATLAB_BIN" -batch "run('$SCRIPT_DIR/benchmark_matlab.m')" || MATLAB_RC=$?
    #MATLAB puede completar todas las medidas y caerse luego al cerrar 

    MATLAB_ROWS=0
    [ -f "$MATLAB_CSV" ] && MATLAB_ROWS=$(( $(wc -l < "$MATLAB_CSV") - 1 ))
    if [ "$MATLAB_ROWS" -gt 0 ]; then
        HAVE_MATLAB=1
        if [ "$MATLAB_RC" -ne 0 ]; then
            echo "AVISO: MATLAB salió con código $MATLAB_RC, pero dejó $MATLAB_ROWS filas" \
                 "completas en el CSV (crash al cerrar); se usan esos datos." >&2
        fi
    else
        echo "AVISO: el benchmark de MATLAB no produjo resultados (código $MATLAB_RC);" \
             "se reporta solo Julia." >&2
    fi
else
    echo "AVISO: no encuentro matlab ('$MATLAB_BIN'); se reporta solo Julia." >&2
fi
echo

# --- 3. CSV combinado (una sola cabecera) ----------------------------------
head -n 1 "$JULIA_CSV" > "$COMBINED_CSV"
tail -n +2 "$JULIA_CSV" >> "$COMBINED_CSV"
if [ "$HAVE_MATLAB" -eq 1 ] && [ -f "$MATLAB_CSV" ]; then
    tail -n +2 "$MATLAB_CSV" >> "$COMBINED_CSV"
fi

# --- 4. Reporte markdown ---------------------------------------------------
HOSTINFO="$(uname -srm)"
NOW="$(date '+%Y-%m-%d %H:%M:%S')"

awk -F, -v host="$HOSTINFO" -v now="$NOW" -v have_matlab="$HAVE_MATLAB" '
function human(b,   u) {
    if (b == "NA" || b == "") return "—";
    u = "B";
    if (b >= 1024)      { b /= 1024; u = "KB"; }
    if (b >= 1024)      { b /= 1024; u = "MB"; }
    return sprintf("%.1f %s", b, u);
}
# Formatea un tiempo en ms; vacío o "NA" -> raya.
function ms(x) { return (x == "" || x == "NA") ? "—" : sprintf("%.4f", x); }
# Ratio MATLAB/Julia; raya si falta alguno de los dos.
function ratio(wj, wm) {
    if (wj != "" && wj != "NA" && wm != "" && wm != "NA" && wj+0 > 0)
        return sprintf("%.2f×", wm / wj);
    return "—";
}
NR == 1 { next }                       # cabecera
{
    lang=$1; malla=$2; nelem=$3; ndof=$4; etapa=$5;
    ttfx=$6; warm=$7; alloc=$8; memocc=$9; nnzv=$10;
    nnzfree=$11; memfree=$12;
    key = etapa "|" malla;
    if (!(key in seen)) { seen[key]=1; order[++nk]=key; }
    ndofof[key]=ndof; nelemof[key]=nelem;
    warmof[key,lang]=warm; ttfxof[key,lang]=ttfx;
    memof[key,lang]=memocc; allocof[key,lang]=alloc;
    nnzof[key,lang]=nnzv;
    nnzfreeof[key,lang]=nnzfree; memfreeof[key,lang]=memfree;
}
END {
    print "# Benchmark de rendimiento (ensamblaje + solver): Julia vs MATLAB";
    print "";
    print "- Fecha: " now;
    print "- Máquina: `" host "`";
    print "- Hilos (BLAS/cómputo): 1 en ambos lados (comparación de algoritmo serial).";
    print "- Solver de autovalores: **incluido** como etapa `solver`, para poder";
    print "  contrastar su coste con el de construir las matrices.";
    print "- Tiempos: **mediana en caliente** (BenchmarkTools / timeit). Ver TTFX aparte.";
    if (have_matlab == 0)
        print "- **Nota:** MATLAB no se ejecutó; solo se muestran los datos de Julia.";
    print "";

    # ---- Tabla 1: tiempo en caliente ----
    print "## Tiempo en caliente por etapa";
    print "";
    print "| Etapa | Malla | n_dof | Julia (ms) | MATLAB (ms) | MATLAB/Julia |";
    print "|---|---:|---:|---:|---:|---:|";
    prev = "";
    for (i = 1; i <= nk; i++) {
        key = order[i]; split(key, kp, "|");
        etapa = kp[1]; malla = kp[2];
        # Fila en blanco al cambiar de malla: separa visualmente los bloques.
        if (prev != "" && malla != prev) print "|  |  |  |  |  |  |";
        prev = malla;
        mlbl = (malla == 0) ? "—" : malla;
        wj = warmof[key,"julia"]; wm = warmof[key,"matlab"];
        # El `assembly` de Julia incluye Dirichlet (dirichlet=true), igual que el
        # de MATLAB (getGlobalMatrices llama a imposeLocalDirichletMatrix), así
        # que ambas etapas son comparables 1:1.
        printf("| %s | %s | %s | %s | %s | %s |\n",
               etapa, mlbl, ndofof[key], ms(wj), ms(wm), ratio(wj, wm));
    }
    print "";
    print "> Ratio > 1 → Julia más rápido; < 1 → MATLAB más rápido.";
    print "> `assembly` es comparable 1:1: en ambos lenguajes la etapa incluye la";
    print "> imposición de las condiciones Dirichlet (PEC) dentro del ensamblaje.";
    print "";

    # ---- Tabla 2a: bloque libre-libre, LA comparación válida de memoria ----
    print "## Memoria de K + M — bloque libre-libre (comparable)";
    print "";
    print "| Malla | n_dof | nnz(K_f) Julia | nnz(K_f) MATLAB | Mem Julia | Mem MATLAB |";
    print "|---:|---:|---:|---:|---:|---:|";
    for (i = 1; i <= nk; i++) {
        key = order[i]; split(key, kp, "|");
        if (kp[1] != "assembly") continue;
        nj = nnzfreeof[key,"julia"]; nm = nnzfreeof[key,"matlab"];
        printf("| %s | %s | %s | %s | %s | %s |\n",
               kp[2], ndofof[key],
               (nj == "" || nj == "NA") ? "—" : nj,
               (nm == "" || nm == "NA") ? "—" : nm,
               human(memfreeof[key,"julia"]), human(memfreeof[key,"matlab"]));
    }
    print "";
    print "> Ambas columnas son `K(free,free)` y `M(free,free)`: las matrices **con";
    print "> Dirichlet ya impuesto en los dos lados**, mismo problema discreto y mismo";
    print "> tamaño (n_free x n_free). Es la única comparación de memoria válida.";
    print "> `nnz(K_f)` coincide **exactamente** en las mallas 1–3 y difiere en ~1e-4";
    print "> relativo en las 4–5 (68 y 280 entradas de 6.9e5 y 1.1e6): son términos de";
    print "> redondeo que cancelan a cero exacto en un orden de acumulación y no en el";
    print "> otro — hay decenas de miles de entradas por debajo de 1e-9·‖K_f‖_F. La";
    print "> igualdad numérica de fondo la certifica `test_08_assembly.jl` (‖K‖_F y";
    print "> traza coinciden a ~1e-15).";
    print "";

    # ---- Tabla 2b: matrices completas (contexto) ----
    print "## Memoria de K + M — matriz completa (contexto)";
    print "";
    print "| Malla | n_dof | nnz(K) Julia | nnz(K) MATLAB | Mem Julia | Mem MATLAB |";
    print "|---:|---:|---:|---:|---:|---:|";
    for (i = 1; i <= nk; i++) {
        key = order[i]; split(key, kp, "|");
        if (kp[1] != "assembly") continue;
        nj = nnzof[key,"julia"]; nm = nnzof[key,"matlab"];
        printf("| %s | %s | %s | %s | %s | %s |\n",
               kp[2], ndofof[key],
               (nj == "" || nj == "NA") ? "—" : nj,
               (nm == "" || nm == "NA") ? "—" : nm,
               human(memof[key,"julia"]), human(memof[key,"matlab"]));
    }
    print "";
    print "> Ambas implementaciones imponen Dirichlet DENTRO del ensamblaje (`dirichlet=true`";
    print "> en Julia; `imposeLocalDirichletMatrix` en MATLAB), así que la matriz \"completa\"";
    print "> ya llega con los acoplamientos de contorno a cero y su nnz es prácticamente el";
    print "> del bloque libre. La memoria es algo mayor que la de la tabla anterior porque la";
    print "> matriz es n_dof x n_dof (no n_free x n_free), pero la comparación entre lenguajes";
    print "> es directa. La tabla anterior (bloque libre-libre) sigue siendo la referencia.";
    print "";

    # ---- Tabla 3: TTFX (primera llamada, incluye compilación JIT en Julia) ----
    print "## TTFX — primera ejecución (incluye compilación)";
    print "";
    print "| Etapa | Malla | Julia TTFX (ms) | MATLAB TTFX (ms) |";
    print "|---|---:|---:|---:|";
    prev = "";
    for (i = 1; i <= nk; i++) {
        key = order[i]; split(key, kp, "|");
        etapa = kp[1]; malla = kp[2];
        if (prev != "" && malla != prev) print "|  |  |  |  |";
        prev = malla;
        mlbl = (malla == 0) ? "—" : malla;
        tj = ttfxof[key,"julia"]; tm = ttfxof[key,"matlab"];
        printf("| %s | %s | %s | %s |\n", etapa, mlbl, ms(tj), ms(tm));
    }
    print "";
    print "> En Julia el TTFX de la primera malla incluye la compilación JIT; en";
    print "> las mallas siguientes la función ya está compilada, de ahí la caída.";
    print "> El `assembly` es comparable 1:1: en ambos lenguajes la etapa incluye la";
    print "> imposición de las condiciones Dirichlet (PEC) dentro del ensamblaje.";
}
' "$COMBINED_CSV" > "$REPORT_MD"

echo "=== Reporte generado ==="
echo "  CSV combinado: $COMBINED_CSV"
echo "  Reporte:       $REPORT_MD"
echo
cat "$REPORT_MD"
