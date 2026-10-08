#!/bin/sh
# Compiles each TikZ figure fig_*.tex in this directory as a standalone document and exports it to PNG (400 dpi).
# Usage: sh build.sh [fig_name ...]      (requires pdflatex with TikZ/pgfplots and pdftoppm)
set -e
cd "$(dirname "$0")"
figs="$*"
[ -z "$figs" ] && figs=$(ls fig_*.tex | sed 's/\.tex$//')
mkdir -p build
for f in $figs; do
  cat > build/$f.tex <<TEX
\documentclass[border=3pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,positioning,calc,decorations.pathreplacing,patterns}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{groupplots,fillbetween}
\setlength{\textwidth}{6.26in}
\begin{document}
\input{../$f}
\end{document}
TEX
  (cd build && pdflatex -interaction=nonstopmode -halt-on-error $f.tex > $f.out 2>&1) || { tail -30 build/$f.log; exit 1; }
  pdftoppm -png -r 400 -singlefile build/$f.pdf $f
  echo "$f.png"
done
