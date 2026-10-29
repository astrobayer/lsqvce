#!/usr/bin/env python3
"""Build the two user documents of lsqvce from their Markdown sources:

  the Manual         dev/docs_src/manual/*.md -> html/lsqvce_manual.html (`doc lsqvce`)
                                                 and lsqvce_manual.pdf
  the API reference  dev/docs_src/api/*.md    -> html/lsqvce_api.html and lsqvce_api.pdf,
                     plus GENERATED chapters so they cannot drift from the code:
                       the help blocks of the user-tier units (fit, Options, Results,
                       report, layout), the engine / kernel / util units from the
                       Contents.m files, and the error-identifier catalogue from
                       dev/docs_src/error_ids.csv.

The theory (lsqvce_theory.pdf) is LaTeX and has its own builder,
dev/docs_src/theory/build.sh.

Run from the toolbox root:  python3 dev/docs_src/build_manual.py
Needs pandoc (html + pdf via xelatex). audit_tree.py check 23 verifies that both
built HTML files equal a fresh build (check mode: build_manual.py check).

THE CONTRACT. Every document is attempted, whatever happened to another. A
document is built ENTIRELY in a temporary directory outside the tree, and its
four outputs -- html/lsqvce_<doc>.html, lsqvce_<doc>.pdf and the tracked
intermediates dev/docs_src/manual/_build/lsqvce_<doc>.md and .html -- replace
the committed ones only when all four built: staged beside their targets (same
file system) and moved with os.replace, the PDF first, because a viewer holding
the PDF open is what makes a target unwritable. Any failure of a document -- a
pandoc or xelatex error, an exception in a generated chapter, a target that
cannot be replaced -- is caught for that document, printed with its name, and
leaves its outputs as they were; the exit code is 1 when any document failed.
Check mode builds the HTML in a temporary directory, compares it with the
committed one and writes nothing in the tree. pandoc's messages are decoded as
UTF-8 with undecodable bytes replaced, whatever the console's code page.
audit_selftest.py drives the build mode with pandoc faked and one fault per case.

Created 2026-09-04 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5.1 -- the docs level of HANDOVER 12.13 (the Manual built from
its Markdown sources; the API reference joined it in WP-124).
Updated 2026-10-17 (CEST, UTC+02:00), WP-146: the contract above -- every
failure caught per document, the four outputs staged and replaced together,
check mode writes nothing, UTF-8 decoding. With assistance from Claude
(Anthropic), model Claude Opus 5.5.
"""
import os, re, shutil, subprocess, sys, csv, io, tempfile

ROOT = os.getcwd()
sys.path.insert(0, os.path.join(ROOT, 'dev', 'docs_src'))
import build_missing_pages as bmp


# a parenthetical that is nothing but a development tag, e.g. "(R5 M1)", "(WP-44)",
# "(R5-D8)": the API states the interface, the record states its history
TAG = re.compile(r' ?\((?:R\d+(?:[- ](?:M|D|WP-?)?\d+\w*)*|WP-\d+\w*)(?: part \d+\w*)?\)')


def help_block(path):
    """The help block of a unit for the API reference: the text up to and including
    its See also paragraph. The dated provenance stamps after it, and bare
    development tags, are history and stay in the source."""
    h1, lines = bmp.read_help(path)
    out = []
    see = False
    for ln in lines:
        if ln.strip().startswith('See also'):
            see = True
        elif see and not ln.strip():
            break
        out.append(TAG.sub('', ln))
    return h1, '\n'.join(out)


def chapter_api():
    out = ['# Help of the user-tier units\n',
           'The help blocks of the six user-tier units, verbatim: `help lsqvce.<name>` shows the same '
           'text and `doc lsqvce.<name>` the same page. Generated from the code at build time.\n']
    for name in ['fit', 'Options', 'Results', 'report', 'layout', 'solveContributions']:
        h1, body = help_block(os.path.join(ROOT, '+lsqvce', name + '.m'))
        out.append('\n## lsqvce.%s\n\n%s\n\n```\n%s\n```\n' % (name, h1.strip(), body.rstrip()))
    return '\n'.join(out)


def chapter_engine():
    out = ['# The engine, kernel and util tiers\n',
           'What `fit` dispatches to and what it is built from: one line per public unit, from the '
           '`Contents.m` files (generated). `doc lsqvce.engine.<name>` opens the full page. Calls go '
           'downward only: the user tier calls the engine, the engine calls kernel and util.\n']
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


def appendix_errors():
    out = ['# Error and warning identifiers\n',
           'Every identifier the package raises (generated from the identifier catalogue, which is kept '
           'closed both ways: nothing is raised that is not listed, nothing is listed that is not raised). '
           'One family is not listed by name: an `lsqvce:makeLoader:<reason>` raised while `fit` builds a '
           'system\'s loader reaches you as `lsqvce:fit:<reason>`, with the original attached as the cause.\n',
           '| identifier | raised by | kind | introduced |\n|---|---|---|---|']
    with open(os.path.join(ROOT, 'dev', 'docs_src', 'error_ids.csv'), encoding='utf-8') as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith('#') or row[0] in ('id', 'identifier'):
                continue
            row = [c.strip().replace('|', '\\|') for c in row] + [''] * 4
            out.append('| `%s` | %s | %s | %s |' % (row[0], row[1].replace(';', ', '), row[2], row[3]))
    return '\n'.join(out) + '\n'



def assemble(src, gens, appendix=None):
    """Concatenate src/*.md in name order; a generated chapter keyed 'NN' goes
    before the first source file whose prefix sorts after it."""
    files = sorted(f for f in os.listdir(src) if f.endswith('.md'))
    parts = []
    seen = set()
    for f in files:
        for g in sorted(gens):
            if g not in seen and g < f[:2]:
                parts.append(gens[g]()); seen.add(g)
        seen.add(f[:2])
        parts.append(open(os.path.join(src, f), encoding='utf-8').read().rstrip('\n') + '\n')
    for g in sorted(gens):
        if g not in seen:
            parts.append(gens[g]()); seen.add(g)
    if appendix:
        parts.append(appendix())
    return '\n\n'.join(parts)


DOCS = {
    'manual': dict(title='User Manual', gens={}, appendix=None),
    'api': dict(title='API Reference', gens={'90': chapter_api, '95': chapter_engine},
                appendix=appendix_errors),
}


class BuildError(Exception):
    """A document's build failed; the message says which step."""


def pandoc(args, what):
    """Run pandoc; its messages decoded as UTF-8 (B6: the locale code page, cp1252
    on Windows, raised on a byte it could not map); a non-zero exit is BuildError."""
    r = subprocess.run(['pandoc'] + args, capture_output=True, encoding='utf-8', errors='replace')
    if r.returncode != 0:
        raise BuildError('%s failed (pandoc exit %d):\n%s' % (what, r.returncode, (r.stderr or '')[-2000:]))


def build_one(name, ver, check):
    """Build one document in a temporary directory. With CHECK, return whether
    the committed HTML differs from the fresh one, writing nothing in the tree.
    Without, replace its four outputs together (the module docstring's
    contract); raise on any failure, before any output was replaced."""
    d = DOCS[name]
    src = os.path.join(ROOT, 'dev', 'docs_src', name)
    stem = 'lsqvce_' + name
    md = ('---\ntitle: "lsqvce %s — %s"\nsubtitle: "least-squares adjustment with variance component estimation"\n'
          '---\n\n' % (ver, d['title'])) + assemble(src, d['gens'], d['appendix'])
    build = os.path.join(ROOT, 'dev', 'docs_src', 'manual', '_build')
    target = os.path.join(ROOT, 'html', stem + '.html')
    with tempfile.TemporaryDirectory() as tmp:          # outside the tree: _build is tracked
        mdpath = os.path.join(tmp, stem + '.md')
        with open(mdpath, 'w', encoding='utf-8') as fh:
            fh.write(md)
        html = os.path.join(tmp, stem + '.html')
        pandoc([mdpath, '-s', '--toc', '--toc-depth=2', '--number-sections', '--metadata', 'lang=en',
                '-c', 'lsqvce_manual.css', '-o', html], 'the HTML build (%s)' % name)
        with open(html, encoding='utf-8') as fh:
            new = fh.read()
        if check:
            cur = ''
            if os.path.exists(target):
                with open(target, encoding='utf-8') as fh:
                    cur = fh.read()
            return cur != new
        pdf = os.path.join(tmp, stem + '.pdf')
        pandoc([mdpath, '--toc', '--toc-depth=2', '--number-sections', '-V', 'geometry:margin=2.2cm',
                '-V', 'fontsize=10pt', '-V', 'colorlinks=true', '--pdf-engine=xelatex',
                '-V', 'fontfamily=mathptmx', '-V', 'mainfont=DejaVu Serif',
                '-V', 'monofont=DejaVu Sans Mono', '-o', pdf], 'the PDF build (%s)' % name)
        # every output built: stage each beside its target, then replace -- the PDF
        # first, the one a viewer locks; a failure leaves every output as it was
        # except when a LATER replace fails, which the caller reports by name
        pairs = [(pdf, os.path.join(ROOT, stem + '.pdf')), (html, target),
                 (mdpath, os.path.join(build, stem + '.md')), (html, os.path.join(build, stem + '.html'))]
        staged = []
        try:
            for built, dst in pairs:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                tmpdst = dst + '.staged'
                shutil.copyfile(built, tmpdst)
                staged.append((tmpdst, dst))
            for k, (tmpdst, dst) in enumerate(staged):
                try:
                    os.replace(tmpdst, dst)
                except OSError as e:
                    done = [os.path.relpath(x, ROOT) for _, x in staged[:k]]
                    err = BuildError('cannot replace %s (%s)%s' % (
                        os.path.relpath(dst, ROOT), e.strerror or e,
                        '' if not done else '; ALREADY replaced: ' + ', '.join(done)))
                    err.partial = bool(done)
                    raise err
        finally:
            for tmpdst, _ in staged:
                if os.path.exists(tmpdst):
                    os.remove(tmpdst)
    print('built html/%s.html and %s.pdf (%s)' % (stem, stem, ver))
    return True


def main(check=False):
    """Build both documents; with check=True return the names whose HTML is stale."""
    ver = re.search(r'^%\s*Version\s+(\S+)', open(os.path.join(ROOT, '+lsqvce', 'Contents.m'),
                    encoding='utf-8').read(), re.M).group(1)
    css = os.path.join(ROOT, 'html', 'lsqvce_manual.css')
    if not check and not os.path.exists(css):
        open(css, 'w', encoding='utf-8').write(CSS)
    if check:
        return [n for n in DOCS if build_one(n, ver, True)]
    written = []
    for n in DOCS:                                       # every one, also after a failure
        try:
            written.append(build_one(n, ver, False))
        except Exception as e:                           # B2: ANY failure is this document's
            left = '' if getattr(e, 'partial', False) else '; its outputs were left as they were'
            print('build failed (%s)%s: %s' % (n, left, e))
            written.append(False)
    return all(written)


CSS = """body{font-family:Georgia,serif;max-width:60em;margin:2em auto;padding:0 1em;line-height:1.45;color:#222}
h1{border-bottom:2px solid #444;padding-bottom:.2em;margin-top:2em}h2{margin-top:1.6em}
code,pre{font-family:Consolas,Menlo,monospace;font-size:.92em}pre{background:#f6f6f6;border:1px solid #ddd;padding:.7em;overflow-x:auto}
table{border-collapse:collapse;margin:1em 0}th,td{border:1px solid #bbb;padding:.3em .6em;vertical-align:top}th{background:#eee}
#TOC{background:#fafafa;border:1px solid #ddd;padding:.5em 1em;margin-bottom:2em}#TOC ul{list-style:none;padding-left:1em}
"""

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'check':
        stale = main(check=True)
        print('stale: %s' % ', '.join(stale) if stale else 'both documents are fresh')
        sys.exit(1 if stale else 0)
    sys.exit(0 if main() else 1)
