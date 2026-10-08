#!/usr/bin/env bash
# build_paper.sh -- build the paper and its source packages into dist/.
#
#   dist/dyadic-blowup.pdf           the compiled paper
#   dist/dyadic-blowup-tex.zip       main.tex with the figures as PNG, and
#                                     their TikZ sources
#   dist/dyadic-blowup-arxiv.tar.gz  main.tex with the figures as PDF
#   dist/dyadic-blowup-physica-d.pdf, dist/dyadic-blowup-physica-d-source.zip
#                                     the same paper in Elsevier's elsarticle
#                                     class for Physica D, assembled from
#                                     main.tex by scripts/build_physica_d.py
#
# The figures are compiled from their TikZ sources by paper/figures/build.sh,
# which also rewrites paper/figures/fig_*.png (with identical content).
# Needs pdflatex with amsart, elsarticle, TikZ and pgfplots; pdftoppm (poppler);
# python3; zip, tar.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAPER="$ROOT/paper"
DIST="$ROOT/dist"
NAME=dyadic-blowup
mkdir -p "$DIST"

latex() { pdflatex -interaction=nonstopmode -halt-on-error "$@" > /dev/null; }
latex3() { latex "$@" && latex "$@" && latex "$@"; }
# check_log LOG [IGNORE]: fail on overfull boxes and undefined or multiply defined
# labels, except lines matching the pattern IGNORE.
check_log() {
  local bad
  bad=$(grep "Overfull\|undefined\|multiply defined" "$1" | grep -v "${2:-^$}" || true)
  if [ -n "$bad" ]; then
    echo "$bad"
    echo "build_paper.sh: fix the warnings above ($1)" >&2
    exit 1
  fi
}
# elsarticle 3.3 itself sets its first-page footer 2.6pt too wide
ELS_FOOTER='Overfull \\hbox (2\.6[0-9]*pt too wide) has occurred while \\output is active'


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

echo "== Physica D version (elsarticle)"
PD="$STAGE/physd/$NAME-physica-d"
mkdir -p "$PD/figures"
python3 "$ROOT/scripts/build_physica_d.py" "$PD/main.tex"
for f in $FIGS; do cp "figures/$f.png" "figures/$f.tex" "$PD/figures/"; done
cp figures/build.sh "$PD/figures/"
cp -r figures/data "$PD/figures/"
( cd "$PD" && latex3 main.tex && check_log main.log "$ELS_FOOTER" \
  && cp main.pdf "$DIST/$NAME-physica-d.pdf" \
  && rm -f main.aux main.log main.out main.pdf main.spl )
rm -f "$DIST/$NAME-physica-d-source.zip"
( cd "$STAGE/physd" && zip -qr "$DIST/$NAME-physica-d-source.zip" "$NAME-physica-d" )

rm -f main.aux main.log main.out main.pdf
ls -l "$DIST"
