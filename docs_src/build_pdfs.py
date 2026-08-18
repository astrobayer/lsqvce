#!/usr/bin/env python3
"""Build the three lsqvce manuals as PDFs, from the SAME authored bodies used for
the HTML (INTERP/REFMAN/USERMAN in build_manuals.py). The bodies are lightweight
HTML; this converts them to LaTeX and compiles with pdflatex. Equations are the
native \\(...\\) / \\[...\\] LaTeX, so the PDF math is real typeset math (not
images). Build-time only: needs pdflatex. No runtime dependency ships.

Run from the toolbox root:  python3 docs_src/build_pdfs.py
Writes lsqvce_user_manual.pdf, lsqvce_reference_manual.pdf,
lsqvce_interpretation_guide.pdf into the toolbox root.
"""
import re, os, subprocess, sys, html as _html

# ---- pull the three authored bodies + their titles from build_manuals.py -----
SRC = open(os.path.join(os.path.dirname(__file__), 'build_manuals.py')).read()

def body_of(name):
    m = re.search(name + r"\s*=\s*r'''(.*?)'''", SRC, re.S)
    if not m:
        raise RuntimeError('body not found: ' + name)
    return m.group(1)

# ---- entity + inline-tag handling -------------------------------------------
ENTITIES = {
    '&amp;': r'\&', '&hellip;': r'\ldots{}', '&ldquo;': '``', '&rdquo;': "''",
    '&lsquo;': '`', '&rsquo;': "'", '&mdash;': '---', '&ndash;': '--',
    '&middot;': r'$\cdot$', '&minus;': r'$-$', '&nbsp;': '~', '&ne;': r'$\ne$',
    '&rarr;': r'$\rightarrow$', '&rsaquo;': r'$\rangle$', '&sect;': r'\S{}',
    '&sigma;': r'$\sigma$', '&sum;': r'$\sum$', '&times;': r'$\times$',
    '&lt;': r'$<$', '&gt;': r'$>$', '&le;': r'$\le$', '&ge;': r'$\ge$',
    '&hArr;': r'$\Leftrightarrow$', '&asymp;': r'$\approx$',
}

# characters that must be escaped in LaTeX text mode are handled by
# esc_text_plain (below), after math/commands have been stashed by inline().

def _esc_verbatimish(c):
    # escape LaTeX specials inside \texttt content; order matters (backslash first)
    c = c.replace('\\', r'\textbackslash{}')
    for ch, rep in [('{', r'\{'), ('}', r'\}'), ('_', r'\_'), ('%', r'\%'),
                    ('#', r'\#'), ('&', r'\&'), ('$', r'\$'),
                    ('^', r'\textasciicircum{}'), ('~', r'\textasciitilde{}')]:
        c = c.replace(ch, rep)
    return c

def inline(s):
    # Convert inline tags to LaTeX, protecting emitted commands + math spans from
    # the text escaper by stashing them as placeholders, escaping the plain text,
    # then restoring.
    stash = []
    def keep(latex):
        stash.append(latex)
        return '\x00%d\x00' % (len(stash) - 1)

    # math spans first (pass through verbatim)
    s = re.sub(r'(\\\(.+?\\\)|\\\[.+?\\\])', lambda m: keep(m.group(1)), s, flags=re.S)
    # code spans
    s = re.sub(r'<code>(.*?)</code>',
               lambda m: keep(r'\texttt{' + _esc_verbatimish(_html.unescape(m.group(1))) + '}'),
               s, flags=re.S)
    # bold / italic (recurse so nested inline works), then stash the result
    s = re.sub(r'<b>(.*?)</b>', lambda m: keep(r'\textbf{' + inline(m.group(1)) + '}'), s, flags=re.S)
    s = re.sub(r'<i>(.*?)</i>', lambda m: keep(r'\emph{' + inline(m.group(1)) + '}'), s, flags=re.S)
    # links -> keep anchor text only (self-contained PDF, no href)
    s = re.sub(r'<a[^>]*>(.*?)</a>', lambda m: inline(m.group(1)), s, flags=re.S)

    s = esc_text_plain(s)                      # escape the remaining plain text
    # restore stashed LaTeX
    s = re.sub(r'\x00(\d+)\x00', lambda m: stash[int(m.group(1))], s)
    return s

def esc_text_plain(s):
    # escape LaTeX specials in text that contains NO math/commands (those were
    # stashed). Entities first, then specials.
    for k, v in ENTITIES.items():
        s = s.replace(k, '\x01%d\x01' % list(ENTITIES).index(k))  # stash entity output
    for ch, rep in [('\\', r'\textbackslash{}'), ('&', r'\&'), ('%', r'\%'),
                    ('#', r'\#'), ('_', r'\_'), ('{', r'\{'), ('}', r'\}'),
                    ('$', r'\$'), ('~', r'\textasciitilde{}'),
                    ('^', r'\textasciicircum{}')]:
        s = s.replace(ch, rep)
    # restore entity replacements
    ent_vals = list(ENTITIES.values())
    s = re.sub(r'\x01(\d+)\x01', lambda m: ent_vals[int(m.group(1))], s)
    return s

# ---- block-level conversion --------------------------------------------------
def convert(body):
    # normalise: strip the class attributes we don't need
    out = []
    # split into top-level blocks by the tags we know
    # process tables specially; everything else line/҂block based
    # Work on a token stream by regex-splitting on block tags.
    # Simpler: iteratively consume recognised blocks.
    i = 0
    text = body
    # handle <pre> ... </pre> as verbatim first (protect content)
    pres = []
    def stash_pre(m):
        pres.append(m.group(1))
        return f'@@PRE{len(pres)-1}@@'
    text = re.sub(r'<pre[^>]*>(.*?)</pre>', stash_pre, text, flags=re.S)

    # tables
    tables = []
    def stash_table(m):
        tables.append(m.group(0))
        return f'@@TABLE{len(tables)-1}@@'
    text = re.sub(r'<table[^>]*>.*?</table>', stash_table, text, flags=re.S)

    lines = []
    # headers
    for raw in re.split(r'(?=<h[123]|<p|<ul|<ol|@@PRE|@@TABLE|<div)', text):
        chunk = raw.strip()
        if not chunk:
            continue
        m = re.match(r'<h1[^>]*>(.*?)</h1>', chunk, re.S)
        if m: lines.append(r'\section*{' + inline(m.group(1)) + '}'); continue
        m = re.match(r'<h2[^>]*>(.*?)</h2>', chunk, re.S)
        if m: lines.append(r'\subsection*{' + inline(m.group(1)) + '}'); continue
        m = re.match(r'<h3[^>]*>(.*?)</h3>', chunk, re.S)
        if m: lines.append(r'\subsubsection*{' + inline(m.group(1)) + '}'); continue
        m = re.match(r'@@PRE(\d+)@@', chunk)
        if m:
            code = _html.unescape(pres[int(m.group(1))])
            lines.append(r'\begin{verbatim}' + '\n' + code.strip('\n') + '\n' + r'\end{verbatim}')
            continue
        m = re.match(r'@@TABLE(\d+)@@', chunk)
        if m:
            lines.append(conv_table(tables[int(m.group(1))]))
            continue
        m = re.match(r'<ul[^>]*>(.*?)</ul>', chunk, re.S)
        if m: lines.append(conv_list(m.group(1), 'itemize')); continue
        m = re.match(r'<ol[^>]*>(.*?)</ol>', chunk, re.S)
        if m: lines.append(conv_list(m.group(1), 'enumerate')); continue
        m = re.match(r'<p[^>]*>(.*?)</p>', chunk, re.S)
        if m:
            # a paragraph may contain a display equation \[...\]
            para = m.group(1)
            para = re.sub(r'\\\[(.+?)\\\]',
                          lambda mm: '\n\\[' + mm.group(1) + '\\]\n', para, flags=re.S)
            lines.append(inline_keep_display(para))
            continue
        # leftover display equation or stray text
        if chunk.startswith(r'\['):
            lines.append(chunk); continue
        # divs / anything else: strip tags, keep text
        chunk = re.sub(r'</?div[^>]*>', '', chunk)
        if chunk.strip():
            lines.append(inline_keep_display(chunk))
    return '\n\n'.join(lines)

def inline_keep_display(s):
    # keep \[...\] as display math, inline the rest
    parts = re.split(r'(\\\[.+?\\\])', s, flags=re.S)
    out = []
    for i, p in enumerate(parts):
        if i % 2 == 1:
            out.append(p)          # display math verbatim
        else:
            out.append(inline(p))
    return ''.join(out)

def conv_list(inner, env):
    items = re.findall(r'<li[^>]*>(.*?)</li>', inner, re.S)
    body = '\n'.join(r'\item ' + inline_keep_display(it.strip()) for it in items)
    return '\\begin{%s}\n%s\n\\end{%s}' % (env, body, env)

def conv_table(tbl):
    rows = re.findall(r'<tr[^>]*>(.*?)</tr>', tbl, re.S)
    grid = []
    for r in rows:
        cells = re.findall(r'<t[hd][^>]*>(.*?)</t[hd]>', r, re.S)
        grid.append([inline_keep_display(c.strip()) for c in cells])
    if not grid:
        return ''
    ncol = max(len(r) for r in grid)
    colspec = 'l' * ncol
    out = [r'\begin{center}\small', r'\begin{tabular}{' + colspec + '}', r'\hline']
    for gi, r in enumerate(grid):
        r = r + [''] * (ncol - len(r))
        out.append(' & '.join(r) + r' \\')
        if gi == 0:
            out.append(r'\hline')
    out += [r'\hline', r'\end{tabular}', r'\end{center}']
    return '\n'.join(out)

# ---- document template -------------------------------------------------------
DOC = r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,bm}
\usepackage[T1]{fontenc}
\usepackage{longtable}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5em}
\title{%(title)s}
\author{lsqvce 1.0 --- consolidated diagnostics + streaming VCE toolbox}
\date{Generated 2026-07-17 with assistance from Claude (Anthropic), Claude Opus 4.8}
\begin{document}
\maketitle
%(body)s
\end{document}
'''

MANUALS = [
    ('USERMAN', 'lsqvce User Manual', 'lsqvce_user_manual'),
    ('REFMAN',  'lsqvce Reference Manual', 'lsqvce_reference_manual'),
    ('INTERP',  'lsqvce Interpretation Guide', 'lsqvce_interpretation_guide'),
]

def main():
    workdir = os.path.abspath('docs_src/_pdfbuild')
    os.makedirs(workdir, exist_ok=True)
    for name, title, stem in MANUALS:
        body = convert(body_of(name))
        tex = DOC % {'title': title, 'body': body}
        texpath = os.path.join(workdir, stem + '.tex')
        open(texpath, 'w').write(tex)
        # two passes for the title/references
        for _ in range(2):
            r = subprocess.run(['pdflatex', '-interaction=nonstopmode',
                                '-halt-on-error', stem + '.tex'],
                               cwd=workdir, capture_output=True, text=True)
        pdf = os.path.join(workdir, stem + '.pdf')
        if not os.path.exists(pdf):
            print('FAILED:', stem)
            print(r.stdout[-2500:])
            sys.exit(1)
        os.replace(pdf, stem + '.pdf')
        print('wrote %s.pdf (%d KB)' % (stem, os.path.getsize(stem + '.pdf') // 1024))

# RETIRED 2026-07-19. The three PDFs this builds -- user manual, reference
# manual, interpretation guide -- were REPLACED by the single master manual
# lsqvce_manual.pdf (docs_src/build_master_pdf.py), and the files themselves were
# deleted from the toolbox root. Running this would recreate three superseded
# documents that no longer agree with the shipped HTML, and it needs a pdflatex
# toolchain that the master manual deliberately does not.
#
# Kept, not deleted, because its HTML->LaTeX conversion is the only record of how
# the old equations were marked up. It refuses to run.
_RETIRED = True

if __name__ == '__main__' and _RETIRED:
    raise SystemExit(
        'build_pdfs.py is RETIRED -- the three PDFs it builds were replaced by '
        'lsqvce_manual.pdf. Build that instead:\n'
        '    python3 docs_src/build_master_pdf.py\n'
        'See HANDOVER 1.3d.')

if __name__ == '__main__':
    main()
