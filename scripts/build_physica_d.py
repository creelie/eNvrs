"""Assemble the Physica D version of the paper as one self-contained .tex file.

    python3 scripts/build_physica_d.py OUT.tex

The text between the markers in paper/main.tex (the shared macros, the main
text and the appendix) is copied unchanged into paper/physica-d/template.tex,
which has Elsevier's elsarticle front matter and declarations. The references
in paper/physica-d/references.tex are put in the order of first citation, as
Elsevier's numbered style asks. The script stops if a cited key has no
reference, if a reference is never cited, or if the two reference lists name
different keys.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper"


def between(text, start, end):
    i, j = text.find(start), text.find(end)
    if i < 0 or j < 0 or j < i:
        sys.exit(f"build_physica_d.py: markers {start!r} ... {end!r} not found in main.tex")
    return text[text.index("\n", i) + 1 : j]


def bibitems(text):
    """Split a list of \\bibitem entries into {key: entry}, keeping their order."""
    parts = re.split(r"(?=\\bibitem\{)", text)
    items = {}
    for p in parts:
        m = re.match(r"\\bibitem\{([^}]+)\}", p)
        if m:
            items[m.group(1)] = p.strip()
    return items


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main_tex = (PAPER / "main.tex").read_text()
    macros = between(main_tex, "%% BEGIN SHARED MACROS", "%% END SHARED MACROS")
    body = between(main_tex, "%% BEGIN BODY", "%% END MAIN TEXT")
    appendix = between(main_tex, "%% END MAIN TEXT", "%% END BODY")

    tex = (PAPER / "physica-d" / "template.tex").read_text()
    for mark, repl in (("%%@MACROS@", macros), ("%%@BODY@", body), ("%%@APPENDIX@", appendix)):
        if tex.count(mark) != 1:
            sys.exit(f"build_physica_d.py: {mark} must occur once in template.tex")
        tex = tex.replace(mark, repl.strip("\n"))

    refs = bibitems((PAPER / "physica-d" / "references.tex").read_text())
    amsart_keys = set(bibitems(main_tex[main_tex.find("\\begin{thebibliography}") :]))
    if set(refs) != amsart_keys:
        sys.exit(f"build_physica_d.py: reference keys differ from main.tex: "
                 f"{sorted(set(refs) ^ amsart_keys)}")

    order = []
    for m in re.finditer(r"\\cite[tp]?\*?(?:\[[^\]]*\])*\{([^}]*)\}", tex):
        for key in (k.strip() for k in m.group(1).split(",")):
            if key not in order:
                order.append(key)
    missing = [k for k in order if k not in refs]
    unused = [k for k in refs if k not in order]
    if missing or unused:
        sys.exit(f"build_physica_d.py: cited without reference {missing}; never cited {unused}")

    if tex.count("%%@REFERENCES@") != 1:
        sys.exit("build_physica_d.py: %%@REFERENCES@ must occur once in template.tex")
    tex = tex.replace("%%@REFERENCES@", "\n\n".join(refs[k] for k in order))
    Path(sys.argv[1]).write_text(tex)


if __name__ == "__main__":
    main()
