#!/usr/bin/env bash
# build_paper.sh -- build the paper and its source packages into dist/.
#
#   dist/dyadic-blowup.pdf           the compiled paper
#   dist/dyadic-blowup-tex.zip       main.tex with the figures as PNG, and
#                                     their TikZ sources
#   dist/dyadic-blowup-arxiv.tar.gz  main.tex with the figures as PDF
#
# The figures are compiled from their TikZ sources by paper/figures/build.sh,
# which also rewrites paper/figures/fig_*.png (with identical content).
# Needs pdflatex with amsart, TikZ and pgfplots; pdftoppm (poppler); zip, tar.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAPER="$ROOT/paper"
DIST="$ROOT/dist"
NAME=dyadic-blowup
mkdir -p "$DIST"

latex() { pdflatex -interaction=nonstopmode -halt-on-error "$@" > /dev/null; }
latex3() { latex "$@" && latex "$@" && latex "$@"; }
check_log() {
  if grep -q "Overfull\|undefined\|multiply defined" "$1"; then
    grep "Overfull\|undefined\|multiply defined" "$1"
    echo "build_paper.sh: fix the warnings above ($1)" >&2
    exit 1
  fi
}

echo "== figures"
sh "$PAPER/figures/build.sh" > /dev/null
FIGS=$(grep -o 'figures/fig_[a-z0-9]*\.png' "$PAPER/main.tex" | sed 's#figures/##; s#\.png$##' | sort -u)

echo "== paper"
cd "$PAPER"
rm -f main.aux main.out
latex3 main.tex
check_log main.log
cp main.pdf "$DIST/$NAME.pdf"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

echo "== tex.zip (figures as PNG)"
mkdir -p "$STAGE/$NAME/figures"
cp main.tex "$STAGE/$NAME/"
for f in $FIGS; do cp "figures/$f.png" "figures/$f.tex" "$STAGE/$NAME/figures/"; done
cp figures/build.sh "$STAGE/$NAME/figures/"
cp -r figures/data "$STAGE/$NAME/figures/"
( cd "$STAGE/$NAME" && latex3 main.tex && check_log main.log \
  && rm -f main.aux main.log main.out main.pdf )
rm -f "$DIST/$NAME-tex.zip"
( cd "$STAGE" && zip -qr "$DIST/$NAME-tex.zip" "$NAME" )

echo "== arXiv tarball (figures as PDF)"
mkdir -p "$STAGE/arxiv/figures"
sed 's#{figures/\(fig_[a-z0-9]*\)\.png}#{figures/\1.pdf}#' main.tex > "$STAGE/arxiv/main.tex"
for f in $FIGS; do cp "figures/build/$f.pdf" "$STAGE/arxiv/figures/"; done
( cd "$STAGE/arxiv" && latex3 main.tex && check_log main.log \
  && rm -f main.aux main.log main.out main.pdf )
tar -czf "$DIST/$NAME-arxiv.tar.gz" -C "$STAGE/arxiv" main.tex figures

rm -f main.aux main.log main.out main.pdf
ls -l "$DIST"
