#!/usr/bin/env bash
# Compila la presentación Beamer.
#
# Se necesitan VARIAS pasadas de pdflatex porque:
#   - El índice (\tableofcontents) se escribe en main.toc durante una pasada
#     y se lee en la siguiente.
#   - La numeración del pie usa \ref{lastmainframe}, que también necesita
#     una pasada previa para resolverse.
# Por eso NO se borran main.toc ni main.aux entre pasadas (solo al final),
# y se hacen 3 pasadas para garantizar que todo quede resuelto.
set -euo pipefail

cd "$(dirname "$0")"

for i in 1 2 3; do
  echo ">> Pasada $i/3"
  pdflatex -interaction=nonstopmode -halt-on-error main.tex > /dev/null
done

# Limpieza de auxiliares (después de compilar, para no romper el índice).
rm -f main.aux main.log main.nav main.out main.snm main.toc main.vrb

echo ">> Listo: $(pwd)/main.pdf"
