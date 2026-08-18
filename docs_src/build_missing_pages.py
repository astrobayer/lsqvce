#!/usr/bin/env python3
"""Generate the per-function HTML doc pages that have no page yet.

WHY THIS EXISTS
    Every lsqvce function is documented by an HTML page under html/, reachable
    from `doc lsqvce.<name>` and listed in helptoc.xml / html/index.html. Twelve
    functions added during the Ellmer / contribution work shipped WITHOUT a page,
    so `doc` failed for them. This script parses the standardised MATLAB help
    block of each such function and emits a page in the existing house style
    (same CSS, same nav/seealso/footer), so the doc set is complete again.

    It is deliberately CONSERVATIVE: it only writes pages that do not exist, and
    it derives everything from the help block, so the page cannot drift from the
    source in ways the help block does not.

USAGE
    python3 docs_src/build_missing_pages.py            # from the toolbox root

Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""

import html
import os
import re
import sys

STYLE = """body { font-family: Helvetica, Arial, sans-serif; font-size: 14px;
       color: #000; background: #fff; margin: 0; padding: 0 20px 40px 20px;
       max-width: 900px; line-height: 1.5; }
h1 { font-size: 24px; font-weight: normal; color: #d9541e;
     border-bottom: 1px solid #d9541e; padding-bottom: 6px; margin-top: 24px; }
h2 { font-size: 18px; font-weight: bold; color: #000; margin-top: 28px; }
h3 { font-size: 15px; font-weight: bold; margin-top: 20px; }
.purpose { font-size: 15px; font-style: italic; color: #444;
           margin: 4px 0 20px 0; }
pre.syntax, pre.example {
  background: #f3f3f3; border: 1px solid #d3d3d3; border-left: 4px solid #d9541e;
  padding: 10px 12px; font-family: Menlo, Consolas, monospace; font-size: 13px;
  overflow-x: auto; white-space: pre-wrap; }
pre.example { border-left-color: #4a90d9; }
code { font-family: Menlo, Consolas, monospace; background: #f3f3f3;
       padding: 1px 4px; border-radius: 2px; }
table { border-collapse: collapse; margin: 12px 0; width: 100%; }
th { background: #eee; text-align: left; padding: 6px 10px;
     border: 1px solid #ccc; font-size: 13px; }
td { padding: 6px 10px; border: 1px solid #ccc; vertical-align: top;
     font-size: 13px; }
td.name { font-family: Menlo, Consolas, monospace; white-space: nowrap;
          width: 1%; }
.warn { background: #fff4e5; border-left: 4px solid #e8912d; padding: 10px 12px;
        margin: 14px 0; }
.danger { background: #fdecea; border-left: 4px solid #c0392b; padding: 10px 12px;
          margin: 14px 0; }
.note { background: #eaf3fb; border-left: 4px solid #4a90d9; padding: 10px 12px;
        margin: 14px 0; }
.nav { font-size: 13px; color: #666; margin: 14px 0 0 0; }
.nav a { color: #0b5ea8; text-decoration: none; }
.nav a:hover { text-decoration: underline; }
a { color: #0b5ea8; }
.seealso { margin-top: 26px; padding-top: 10px; border-top: 1px solid #ddd;
           font-size: 13px; }
.footer { margin-top: 34px; padding-top: 10px; border-top: 1px solid #ddd;
          font-size: 12px; color: #888; }"""

FOOTER = ("lsqvce 1.0 &mdash; consolidated diagnostics + streaming VCE toolbox. "
          "Page generated from the function help block with assistance from "
          "Claude (Anthropic), Claude Opus 4.8, 2026-07-18.")

# Section headings recognised inside a help block. Everything that is not one of
# these is treated as description prose.
SECTION = re.compile(
    r'^(SYNTAX|INPUTS?|OUTPUTS?|OPTIONS?|NAME-VALUE.*|ERRORS?|EXAMPLES?|NOTES?|'
    r'REFERENCES?|ALGORITHM|CONVENTIONS?|WHY.*|WHAT.*|COST.*|SEE ALSO|'
    r'BEST PRACTICES?|ACCURACY|LIMITATIONS?|CONTRACT.*|THE .*)$')


def read_help(path):
    """Return (h1_line, [help lines]) for a .m file, stripping the leading '%'."""
    lines = []
    started = False
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            s = raw.rstrip('\n')
            if not started:
                # skip until the function/classdef declaration
                if re.match(r'^\s*(function|classdef)\b', s):
                    started = True
                continue
            if s.strip().startswith('%'):
                lines.append(re.sub(r'^\s*%\s?', '', s))
            elif s.strip() == '':
                if lines:
                    break
            else:
                break
    if not lines:
        return '', []
    return lines[0], lines


def parse(lines):
    """Split a help block into an ordered list of (heading|None, [lines])."""
    blocks = []
    cur_head, cur = None, []
    for ln in lines[1:]:                     # lines[0] is the H1 summary
        stripped = ln.strip()
        # A heading is an ALL-CAPS-ish line at low indent with no trailing text.
        if stripped and not ln.startswith('    ') and SECTION.match(stripped):
            if cur:
                blocks.append((cur_head, cur))
            cur_head, cur = stripped, []
        else:
            cur.append(ln)
    if cur:
        blocks.append((cur_head, cur))
    return blocks


def trim(block):
    """Drop leading/trailing blank lines."""
    while block and not block[0].strip():
        block.pop(0)
    while block and not block[-1].strip():
        block.pop()
    return block


def dedent(block):
    """Remove the common leading indent so <pre> content lines up."""
    ind = [len(l) - len(l.lstrip()) for l in block if l.strip()]
    k = min(ind) if ind else 0
    return [l[k:] if len(l) >= k else l for l in block]


def linkify_seealso(text):
    """Turn 'LSQVCE.FOO, LSQVCE.BAR' into links to the corresponding pages."""
    out = []
    for tok in re.split(r',\s*', text.strip().rstrip('.')):
        t = tok.strip()
        if not t:
            continue
        m = re.match(r'^(LSQVCE\.[A-Za-z_][\w.]*)$', t, re.I)
        if m:
            name = m.group(1)
            page = 'lsqvce.' + name.split('.', 1)[1]
            out.append('<a href="%s.html">%s</a>' % (html.escape(page),
                                                     html.escape(page)))
        else:
            out.append('<code>%s</code>' % html.escape(t))
    return ', '.join(out)


def render(name, h1, blocks):
    """Build the full HTML page for one function."""
    title = 'lsqvce.' + name
    purpose = re.sub(r'^\s*[A-Z][A-Z0-9_.]*\s+', '', h1).strip()
    if not purpose:
        purpose = h1.strip()

    body = []
    seealso_html = ''
    for head, raw in blocks:
        blk = dedent(trim(list(raw)))
        if not blk:
            continue
        joined = '\n'.join(blk)
        if head and head.upper().startswith('SEE ALSO'):
            seealso_html = linkify_seealso(joined)
            continue
        if head is None:
            # free description prose before the first heading
            body.append('<h2>Description</h2>')
            body.append('<pre class="example">%s</pre>' % html.escape(joined))
        else:
            cls = 'syntax' if head.upper().startswith('SYNTAX') else 'example'
            body.append('<h3>%s</h3>' % html.escape(head.title()
                                                   if head.isupper() else head))
            body.append('<pre class="%s">%s</pre>' % (cls, html.escape(joined)))

    if not seealso_html:
        seealso_html = '<a href="index.html">lsqvce function index</a>'

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>{title}</title><style>
{style}
</style></head>
<body>
<div class="nav"><a href="index.html">lsqvce</a> &rsaquo; {title}</div>
<h1>{title}</h1>
<div class="purpose">{purpose}</div>
{body}
<div class="seealso"><b>See also:</b> {seealso}</div>
<div class="footer">{footer}</div>
</body></html>
""".format(title=html.escape(title), style=STYLE,
           purpose=html.escape(purpose), body='\n'.join(body),
           seealso=seealso_html, footer=FOOTER)


def main():
    root = os.getcwd()
    pkg = os.path.join(root, '+lsqvce')
    if not os.path.isdir(pkg):
        sys.exit('run from the toolbox root (the folder holding +lsqvce)')

    written = []
    targets = []
    for fn in sorted(os.listdir(pkg)):
        if fn.endswith('.m') and fn != 'Contents.m':
            targets.append((fn[:-2], os.path.join(pkg, fn),
                            'lsqvce.%s.html' % fn[:-2]))
    # +util is a sub-package and its pages are named lsqvce.util.<name>.html.
    # It was originally missed here, which is how util.verbosity shipped without
    # a page; scan it explicitly rather than relying on the flat listing.
    upkg = os.path.join(pkg, '+util')
    if os.path.isdir(upkg):
        for fn in sorted(os.listdir(upkg)):
            if fn.endswith('.m') and fn != 'Contents.m':
                targets.append(('util.' + fn[:-2], os.path.join(upkg, fn),
                                'lsqvce.util.%s.html' % fn[:-2]))

    for name, src, pagename in targets:
        page = os.path.join(root, 'html', pagename)
        if os.path.exists(page):
            continue
        h1, lines = read_help(src)
        if not lines:
            print('  no help block, skipped: %s' % name)
            continue
        with open(page, 'w', encoding='utf-8') as fh:
            fh.write(render(name, h1, parse(lines)))
        written.append(name)
        print('  wrote html/%s' % pagename)

    print('\n%d page(s) written.' % len(written))
    return written


if __name__ == '__main__':
    main()
