#!/usr/bin/env bash
#
# build.sh — build paper/main.tex and paper/main.pdf from paper/main.md
#            and refresh MANIFEST.sha256.
#
# Source of truth:
#     paper/references.bib    — bibliography (edit this for any reference change)
#     paper/main.md           — manuscript body (everything before "## References")
#
# Pipeline:
#   1. Regenerate the "## References" section of main.md from references.bib
#      (the bib is the single source of truth; the Markdown reference list
#      below "## References" is overwritten on every build).
#   2. Pandoc converts main.md -> main.tex.
#   3. Small sed/python patches apply LaTeX-only fixes that the Markdown
#      source cannot carry (resizebox for the §1.1 underbrace chain;
#      \path{...} + \small for the §8.2 script-to-theorem table).
#   4. pdflatex is run twice so that cross-references and table widths
#      settle.
#   5. MANIFEST.sha256 is refreshed over every user-visible artefact.
#
# Run from the repository root:
#     ./build.sh

set -euo pipefail

cd "$(dirname "$0")"      # repo root
cd paper                  # pandoc/pdflatex operate inside paper/

echo "[1/6] regenerate ## References in main.md from references.bib"
python3 - <<'PY'
import re
import sys
from pathlib import Path

BIB = Path("references.bib")
MD = Path("main.md")


def parse_bib(text):
    """Return list of (key, entry_type, fields) preserving order."""
    entries = []
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)\n\}\s*",
                         text, re.DOTALL):
        etype = m.group(1).lower()
        key = m.group(2)
        body = m.group(3)
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"[^\"]*\")\s*,?",
                              body, re.DOTALL):
            name = fm.group(1).lower()
            val = fm.group(2).strip()
            if val.startswith("{") and val.endswith("}"):
                val = val[1:-1]
            elif val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            val = re.sub(r"\s+", " ", val).strip()
            # Strip BibTeX case-protection braces: {J}acobi -> Jacobi
            val = re.sub(r"\{([^{}]*)\}", r"\1", val)
            # Convert common LaTeX accents to Unicode.
            val = (val
                   .replace('\\"a', 'ä').replace('\\"o', 'ö').replace('\\"u', 'ü')
                   .replace("\\'a", 'á').replace("\\'e", 'é').replace("\\'i", 'í')
                   .replace("\\'o", 'ó').replace("\\'u", 'ú')
                   .replace('\\"A', 'Ä').replace('\\"O', 'Ö').replace('\\"U', 'Ü'))
            fields[name] = val
        entries.append((key, etype, fields))
    return entries


def format_authors(s):
    """BibTeX 'Last1, F1 and Last2, F2 and ...' -> 'F1 Last1, F2 Last2 and FN LastN'."""
    parts = [p.strip() for p in s.split(" and ")]
    out = []
    for p in parts:
        if "," in p:
            last, first = [x.strip() for x in p.split(",", 1)]
            out.append(f"{first} {last}")
        else:
            out.append(p)
    if len(out) == 1:
        return out[0]
    if len(out) == 2:
        return f"{out[0]} and {out[1]}"
    return ", ".join(out[:-1]) + " and " + out[-1]


def format_entry(key, etype, f):
    """Render one Markdown reference line: [KEY] Authors, *Title*, venue, year."""
    authors = format_authors(f.get("author", ""))
    title = f.get("title", "").replace("\n", " ")
    parts = [f"[{key}] {authors}, *{title}*"]

    if etype == "article":
        journal = f.get("journal", "")
        vol = f.get("volume", "")
        pages = f.get("pages", "").replace("--", "–")
        year = f.get("year", "")
        venue = journal
        if vol:
            venue += f" **{vol}**"
        if year:
            venue += f" ({year})"
        if pages:
            venue += f", {pages}"
        parts.append(venue)
    elif etype == "incollection":
        booktitle = f.get("booktitle", "")
        pub = f.get("publisher", "")
        addr = f.get("address", "")
        year = f.get("year", "")
        pages = f.get("pages", "").replace("--", "–")
        venue = f"in: *{booktitle}*"
        tail = ", ".join(x for x in (pub, addr, year) if x)
        if tail:
            venue += f", {tail}"
        if pages:
            venue += f", {pages}"
        parts.append(venue)
    elif etype == "book":
        series = f.get("series", "")
        vol = f.get("volume", "")
        pub = f.get("publisher", "")
        addr = f.get("address", "")
        year = f.get("year", "")
        edition = f.get("edition", "")
        venue_bits = []
        if series:
            seg = series
            if vol:
                seg += f" **{vol}**"
            venue_bits.append(seg)
        if edition:
            suffix = {"1": "st", "2": "nd", "3": "rd"}.get(edition, "th")
            venue_bits.append(f"{edition}{suffix} ed.")
        if pub:
            venue_bits.append(pub)
        if addr:
            venue_bits.append(addr)
        if year:
            venue_bits.append(year)
        parts.append(", ".join(venue_bits))
    else:
        # Fallback: dump remaining fields as a flat list.
        rest = [v for k, v in f.items() if k not in {"author", "title"}]
        parts.append(", ".join(rest))

    line = ", ".join(parts) + "."
    return line


bib_text = BIB.read_text()
entries = parse_bib(bib_text)
if not entries:
    sys.exit("ERROR: no entries parsed from references.bib")

generated = "## References\n\n"
for key, etype, fields in entries:
    generated += format_entry(key, etype, fields) + "\n\n"

md_text = MD.read_text()
m_start = re.search(r"^## References\b", md_text, re.MULTILINE)
if not m_start:
    sys.exit("ERROR: '## References' heading not found in main.md")
# Find the next top-level heading after ## References (or EOF if none).
m_next = re.search(r"^## ", md_text[m_start.end():], re.MULTILINE)
if m_next:
    end_idx = m_start.end() + m_next.start()
else:
    end_idx = len(md_text)
new_md = md_text[:m_start.start()] + generated + md_text[end_idx:]
MD.write_text(new_md)
print(f"  wrote {len(entries)} entries to main.md")
PY

echo "[2/6] pandoc main.md -> main.tex"
pandoc main.md -o main.tex --standalone --bibliography=references.bib

# ---------------------------------------------------------------------
# Patch 1: Pandoc renders the Markdown [text] in math contexts as
#          \texttt{ ... } and wraps brackets as {[}...{]}. Normalise both.
# ---------------------------------------------------------------------
echo "[3/6] applying LaTeX-only patches"
sed -i 's/{\[}/[/g; s/{\]}/]/g' main.tex

# ---------------------------------------------------------------------
# Patch 2: §1.1 chain — wrap the wide underbrace cascade in \resizebox
#          so it fits the text width.
# ---------------------------------------------------------------------
python3 - <<'PY'
import re
path = "main.tex"
src = open(path).read()
pat = re.compile(
    r"(\\\[)\s*\n"
    r"(\\underbrace\{K\(5,2\)\}_\{\\text\{Petersen\}\}\s*\n"
    r"(?:.*\n)*?"
    r"\\underbrace\{\\text\{icosahedral config\}\}_\{\\text\{rigid\}\})\s*\n"
    r"(\\\])",
)
def wrap(m):
    return f"{m.group(1)}\n\\resizebox{{\\linewidth}}{{!}}{{$\n{m.group(2)}\n$}}\n{m.group(3)}"
new = pat.sub(wrap, src, count=1)
if new == src:
    print("  (note: §1.1 resizebox pattern not matched -- skipped)")
else:
    open(path, "w").write(new)
PY

# ---------------------------------------------------------------------
# Patch 3: §8.2 script-to-theorem mapping table.
#          - Wrap in {\small ... } so the cell text is one size smaller.
#          - Replace \texttt{xxx\_yyy.py} with \path{xxx_yyy.py} so the
#            file names break at underscores (\texttt is non-breakable).
#          - Restore one EXACT_RANK_MODE textt cell that we want as a
#            literal inline keyword, not as a path.
# ---------------------------------------------------------------------
python3 - <<'PY'
import re
path = "main.tex"
src = open(path).read()

# Locate the §8.2 table block (heading -> end of longtable).
heading_marker = "8.2 Script"
end_marker = "\\end{longtable}"
i = src.find(heading_marker)
if i < 0:
    print("  (note: §8.2 heading not found -- skipped)")
    raise SystemExit
j = src.find(end_marker, i)
if j < 0:
    print("  (note: §8.2 longtable end not found -- skipped)")
    raise SystemExit
j_end = j + len(end_marker)

block = src[i:j_end]

# 1) Insert \small before \begin{longtable} (and a matching } after).
block = block.replace("\\begin{longtable}", "{\\small\n\\begin{longtable}", 1)
block = block + "\n}"

# 2) Replace \texttt{filename.py} with \path{filename.py} (also
#    unescaping \_ to _ inside that argument).
def texttt_to_path(m):
    inside = m.group(1).replace("\\_", "_")
    return r"\path{" + inside + "}"
block = re.sub(r"\\texttt\{([^{}]*\.py)\}", texttt_to_path, block)

# 3) Re-insert the EXACT_RANK_MODE inline literal as plain \texttt
#    (was accidentally caught by the \path{...py} pass? No, EXACT_RANK_MODE
#    has no .py; just guard for completeness).
block = block.replace(r"\path{EXACT_RANK_MODE_=_True}",
                      r"\texttt{EXACT\_RANK\_MODE = True}")

src = src[:i] + block + src[j_end:]
open(path, "w").write(src)
PY

# ---------------------------------------------------------------------
# Patch 4: stray Unicode glyphs that pdflatex cannot handle.
# ---------------------------------------------------------------------
sed -i 's/120°/120^\\circ/g' main.tex   # no-op if already done in main.md

# ---------------------------------------------------------------------
# Patch 5: convert in-text [KEY] / [K1, K2, ...] citations to \cite{...}
#          and replace the Pandoc-rendered "References" subsection with
#          a proper \bibliography{references} call. This wires up the
#          standard LaTeX/BibTeX pipeline so the linter sees every .bib
#          entry being used, and the bibliography is generated from the
#          .bib file directly.
# ---------------------------------------------------------------------
python3 - <<'PY'
import re

path = "main.tex"
src = open(path).read()

# Parse the bib once to get the canonical key list.
bib = open("references.bib").read()
keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", bib))


def cite_sub(m):
    inside = m.group(1)
    parts = [p.strip() for p in inside.split(",")]
    if parts and all(p in keys for p in parts):
        return "\\cite{" + ",".join(parts) + "}"
    return m.group(0)  # leave untouched if not a list of bib keys


# Match [KEY] or [K1, K2, ...]; allow line breaks inside the bracket.
src = re.sub(
    r"\[\s*((?:[A-Za-z][A-Za-z0-9]*)(?:\s*,\s*[A-Za-z][A-Za-z0-9]*)*)\s*\]",
    cite_sub,
    src,
    flags=re.DOTALL,
)

# Replace the Pandoc-rendered References subsection (everything from
# "\subsection{References}" up to "\end{document}") with the standard
# bibliography call.
pat = re.compile(
    r"\\subsection\{References\}\\label\{references\}.*?(?=\\end\{document\})",
    re.DOTALL,
)
src = pat.sub(
    "\\\\bibliographystyle{plain}\n\\\\bibliography{references}\n\n",
    src,
)

open(path, "w").write(src)
PY

# ---------------------------------------------------------------------
# Patch 6: swap raster figure paths in \includegraphics to the vector
#          PDF versions. main.md keeps .png so figures still render on
#          GitHub's Markdown view; the LaTeX/arXiv build uses .pdf for
#          true vector quality (matplotlib output is vector-native).
# ---------------------------------------------------------------------
sed -i -E 's|(\\includegraphics(\[[^]]*\])?\{figures/[A-Za-z0-9_.-]+)\.png\}|\1.pdf}|g' main.tex

# ---------------------------------------------------------------------
# Compile: pdflatex -> bibtex -> pdflatex -> pdflatex
# (extra passes so longtable widths, cross-references, and the
# bibliography settle.)
# ---------------------------------------------------------------------
echo "[4/6] pdflatex (pass 1) + bibtex + pdflatex (pass 2 + 3)"
pdflatex -interaction=nonstopmode main.tex >/tmp/pdflatex1.log
bibtex main >/tmp/bibtex.log
pdflatex -interaction=nonstopmode main.tex >/tmp/pdflatex2.log
echo "[5/6] pdflatex (final pass)"
pdflatex -interaction=nonstopmode main.tex >/tmp/pdflatex3.log

# Report any overfull boxes (best-effort).
overfull=$(grep -c "Overfull" main.log || true)
if [ "${overfull:-0}" -gt 0 ]; then
    echo "WARNING: ${overfull} overfull box(es) in main.log; check before submission."
else
    echo "OK: no overfull boxes."
fi

# ---------------------------------------------------------------------
# Step 6: refresh the repository-root MANIFEST.sha256 over every
#         user-visible artefact (sources, generated tex/pdf, figures,
#         scripts, supplement logs, viz, top-level metadata).
# ---------------------------------------------------------------------
echo "[6/6] refresh MANIFEST.sha256"
(
    cd ..
    (
        find paper scripts supplement viz -type f \( \
            -name "*.py" -o -name "*.md" -o -name "*.tex" -o -name "*.pdf" \
            -o -name "*.png" -o -name "*.html" -o -name "*.log" \
            -o -name "*.txt" -o -name "*.cff" -o -name "*.yml" \
            -o -name "*.bib" -o -name "*.sh" \) | sort
        printf '%s\n' README.md requirements.txt LICENSE LICENSE-CC-BY-4.0 \
                      CITATION.cff build.sh
    ) | xargs sha256sum > MANIFEST.sha256
    echo "  $(wc -l < MANIFEST.sha256) hashes"
)

echo
ls -la main.pdf
