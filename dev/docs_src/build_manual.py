#!/usr/bin/env python3
"""Build the lsqvce Manual -- the single documentation authority -- from
dev/docs_src/manual/*.md plus three GENERATED chapters, into
html/lsqvce_manual.html (the page `doc lsqvce` opens) and lsqvce_manual.pdf.

Generated (so they cannot drift from the code):
  05  the user-tier API: the full help block of fit, Options, Results, report, layout
  10  the engine / kernel / util tiers: one line per unit from the Contents.m files
  12  the migration table v4 -> v5, from RELEASE_NOTES.md
  A   the error-identifier catalogue, from dev/docs_src/error_ids.csv

Run from the toolbox root:  python3 dev/docs_src/build_manual.py
Needs pandoc (html + pdf via pdflatex). audit_tree.py check 23 verifies the built
HTML equals a fresh build.
"""
import os, re, subprocess, sys, csv, io

ROOT = os.getcwd()
SRC = os.path.join(ROOT, 'dev', 'docs_src', 'manual')
sys.path.insert(0, os.path.join(ROOT, 'dev', 'docs_src'))
import build_missing_pages as bmp


def help_block(path):
    h1, lines = bmp.read_help(path)
    return h1, '\n'.join(lines)


def chapter_api():
    out = ['# The user tier: fit, Options, Results, report, layout\n',
           'The help blocks of the five user-tier units, verbatim (`help lsqvce.<name>` shows the same text; '
           '`doc lsqvce.<name>` the same page). This chapter is generated from the code at build time.\n']
    for name in ['fit', 'Options', 'Results', 'report', 'layout']:
        h1, body = help_block(os.path.join(ROOT, '+lsqvce', name + '.m'))
        out.append('\n## lsqvce.%s\n\n%s\n\n```\n%s\n```\n' % (name, h1.strip(), body.rstrip()))
    return '\n'.join(out)


def chapter_engine():
    out = ['# The engine, kernel and util tiers\n',
           'What `fit` dispatches to and what it is built from. One line per public unit, from the '
           '`Contents.m` files (generated). `doc lsqvce.engine.<name>` opens the full page; the tier rule '
           '(HANDOVER, tiers.csv): the user tier calls engine, engine calls kernel and util, never upward.\n']
    for pkg, title in [('+engine', 'lsqvce.engine — drivers, accumulator, loaders, sweeps, diagnostics'),
                       ('+kernel', 'lsqvce.kernel — the numerical kernels'),
                       ('+util', 'lsqvce.util — plumbing')]:
        p = os.path.join(ROOT, '+lsqvce', pkg, 'Contents.m')
        if not os.path.exists(p):
            continue
        out.append('\n## %s\n\n| unit | what |\n|---|---|' % title)
        for line in open(p, encoding='utf-8'):
            m = re.match(r'%\s{2,}(\w+)\s+-\s+(.*)', line)
            if m:
                out.append('| `%s` | %s |' % (m.group(1), m.group(2).strip()))
    return '\n'.join(out) + '\n'


def chapter_migration():
    t = open(os.path.join(ROOT, 'RELEASE_NOTES.md'), encoding='utf-8').read()
    i = t.find('### Migration table (v4 → v5)')
    j = t.find('\n## ', i)
    tbl = t[i:j if j > 0 else None]
    tbl = tbl.replace('### Migration table (v4 → v5)', '## The migration table (from RELEASE_NOTES.md)')
    head = ('# Migration from 4.x to 5\n\nEvery 4.x declaration still fits; three things were added '
            '(priors, per-stream policy, `OnPass`), the option space was cut to 14 knobs + 5 groups, '
            'and every retired name has the successor below. Numbers may differ from 4.0 because the '
            'separate-mode estimator (an AR law fitted to a biased raw-residual ACF) is gone; the '
            'acceptance of its successor is the physical oracle, not equality to 4.0.\n\n')
    return head + tbl + '\n'


def appendix_errors():
    out = ['# Appendix A — error and warning identifiers\n',
           'Every identifier the package raises, from `dev/docs_src/error_ids.csv` (the catalogue the audit '
           'keeps closed both ways: nothing raised that is not listed, nothing listed that is not raised).\n',
           '| identifier | where / what |\n|---|---|']
    with open(os.path.join(ROOT, 'dev', 'docs_src', 'error_ids.csv'), encoding='utf-8') as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith('#') or row[0] in ('id', 'identifier'):
                continue
            out.append('| `%s` | %s |' % (row[0], ' — '.join(c.strip() for c in row[1:] if c.strip()).replace('|', '\\|')))
    return '\n'.join(out) + '\n'


def assemble():
    files = sorted(f for f in os.listdir(SRC) if f.endswith('.md'))
    parts = []
    gens = {'05': chapter_api, '10': chapter_engine, '12': chapter_migration}
    seen = set()
    for f in files:
        for g in sorted(gens):                       # generated chapters go BEFORE the first later file
            if g not in seen and g < f[:2]:
                parts.append(gens[g]()); seen.add(g)
        seen.add(f[:2])
        parts.append(open(os.path.join(SRC, f), encoding='utf-8').read().rstrip('\n') + '\n')
    for g in sorted(gens):
        if g not in seen:
            parts.append(gens[g]()); seen.add(g)
    parts.append(appendix_errors())
    return '\n\n'.join(parts)


def main(check=False):
    ver = re.search(r'^%\s*Version\s+(\S+)', open(os.path.join(ROOT, '+lsqvce', 'Contents.m'), encoding='utf-8').read(), re.M).group(1)
    md = ('---\ntitle: "lsqvce %s — User Manual"\nsubtitle: "least-squares adjustment with variance component estimation"\n'
          'date: "built from the sources; the single documentation authority"\n---\n\n' % ver) + assemble()
    build = os.path.join(ROOT, 'dev', 'docs_src', 'manual', '_build')
    os.makedirs(build, exist_ok=True)
    mdpath = os.path.join(build, 'lsqvce_manual.md')
    open(mdpath, 'w', encoding='utf-8').write(md)
    html = os.path.join(build, 'lsqvce_manual.html')
    subprocess.run(['pandoc', mdpath, '-s', '--toc', '--toc-depth=2', '--number-sections',
                    '--metadata', 'lang=en', '-c', 'lsqvce_manual.css', '-o', html], check=True)
    css = os.path.join(ROOT, 'html', 'lsqvce_manual.css')
    if not os.path.exists(css):
        open(css, 'w', encoding='utf-8').write(CSS)
    target = os.path.join(ROOT, 'html', 'lsqvce_manual.html')
    new = open(html, encoding='utf-8').read()
    if check:
        cur = open(target, encoding='utf-8').read() if os.path.exists(target) else ''
        return cur != new
    open(target, 'w', encoding='utf-8').write(new)
    pdf = os.path.join(ROOT, 'lsqvce_manual.pdf')
    r = subprocess.run(['pandoc', mdpath, '--toc', '--toc-depth=2', '--number-sections',
                        '-V', 'geometry:margin=2.2cm', '-V', 'fontsize=10pt', '-V', 'colorlinks=true',
                        '--pdf-engine=xelatex', '-V', 'fontfamily=mathptmx', '-V', 'mainfont=DejaVu Serif', '-V', 'monofont=DejaVu Sans Mono', '-o', pdf], capture_output=True, text=True)
    if r.returncode != 0:
        print('PDF build failed:\n' + r.stderr[-2000:])
        return False
    print('built html/lsqvce_manual.html and lsqvce_manual.pdf (%s)' % ver)
    return True


CSS = """body{font-family:Georgia,serif;max-width:60em;margin:2em auto;padding:0 1em;line-height:1.45;color:#222}
h1{border-bottom:2px solid #444;padding-bottom:.2em;margin-top:2em}h2{margin-top:1.6em}
code,pre{font-family:Consolas,Menlo,monospace;font-size:.92em}pre{background:#f6f6f6;border:1px solid #ddd;padding:.7em;overflow-x:auto}
table{border-collapse:collapse;margin:1em 0}th,td{border:1px solid #bbb;padding:.3em .6em;vertical-align:top}th{background:#eee}
#TOC{background:#fafafa;border:1px solid #ddd;padding:.5em 1em;margin-bottom:2em}#TOC ul{list-style:none;padding-left:1em}
"""

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'check':
        sys.exit(1 if main(check=True) else 0)
    main()
