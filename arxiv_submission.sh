#!/usr/bin/env bash
#
# arxiv_submission.sh — bundle the TeX source for arXiv upload.
#
# Run AFTER ./build.sh so paper/main.tex is fresh. Produces
# arxiv_submission.tar.gz at the repo root, containing exactly the
# files arXiv needs:
#
#     main.tex
#     references.bib
#     figures/<files referenced in main.tex>
#
# Auxiliary files (.aux, .log, .bbl, .blg, .out, .toc) and the
# pre-built main.pdf are excluded — arXiv recompiles from source.
#
# Before archiving, the staged tree is compiled standalone
# (pdflatex -> bibtex -> pdflatex -> pdflatex) to confirm that the
# bundle is self-contained and that arXiv will be able to build it.
#
# Run from the repository root:
#     ./arxiv_submission.sh

set -euo pipefail

cd "$(dirname "$0")"

OUT="arxiv_submission.tar.gz"
STAGE=$(mktemp -d -t arxiv-XXXXXX)
trap 'rm -rf "$STAGE"' EXIT

if [ ! -f paper/main.tex ]; then
    echo "ERROR: paper/main.tex not found. Run ./build.sh first." >&2
    exit 1
fi

echo "[1/4] staging files in $STAGE"
cp paper/main.tex "$STAGE/"
cp paper/references.bib "$STAGE/"
mkdir "$STAGE/figures"

# Pull only the figure files actually referenced in main.tex.
mapfile -t FIGS < <(grep -oE 'figures/[A-Za-z0-9_.-]+\.(png|pdf)' paper/main.tex | sort -u)
if [ "${#FIGS[@]}" -eq 0 ]; then
    echo "ERROR: no figures referenced in main.tex" >&2
    exit 1
fi
for f in "${FIGS[@]}"; do
    if [ ! -f "paper/$f" ]; then
        echo "ERROR: paper/main.tex references $f but the file is missing." >&2
        exit 1
    fi
    cp "paper/$f" "$STAGE/$f"
    echo "  + $f"
done

echo "[2/4] verifying the staged tree compiles standalone"
(
    cd "$STAGE"
    pdflatex -interaction=nonstopmode main.tex >/tmp/arxiv_pdflatex1.log
    bibtex main >/tmp/arxiv_bibtex.log
    pdflatex -interaction=nonstopmode main.tex >/tmp/arxiv_pdflatex2.log
    pdflatex -interaction=nonstopmode main.tex >/tmp/arxiv_pdflatex3.log
)
overfull=$(grep -c "Overfull" "$STAGE/main.log" || true)
if [ "${overfull:-0}" -gt 0 ]; then
    echo "  WARNING: ${overfull} overfull box(es) in the staged compile."
else
    echo "  OK: standalone compile clean."
fi

echo "[3/4] stripping auxiliary files before archiving"
rm -f "$STAGE"/main.aux \
      "$STAGE"/main.log \
      "$STAGE"/main.bbl \
      "$STAGE"/main.blg \
      "$STAGE"/main.out \
      "$STAGE"/main.toc \
      "$STAGE"/main.pdf

echo "[4/4] creating $OUT"
tar -czf "$OUT" -C "$STAGE" .

echo
ls -la "$OUT"
echo
echo "Contents:"
tar -tzf "$OUT" | sort
