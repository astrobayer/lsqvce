#!/usr/bin/env python3
"""refactor_fluent.py -- R4-D2 / WP-64b: rewrite every fluent configurator chain
of lsqvce.engine.adjustVCE into constructor Name-Value option groups.

    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20) ...
         .varianceComponents('MaxIter', 4);
becomes, on the statement that BUILT P,
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 2, 'GroupBy', 'obsType', ...
                  'Shape', struct('MaxOrder', 20), 'Scale', struct('MaxIter', 4));
and the chain statement disappears (or keeps its non-configurator tail,
e.g. `.run()`, `.empiricalStochastic(...)`). Inline constructions get the
groups inserted into the constructor call. `.joint(...)` also sets
'VceMode', 'joint' (the deleted method did).

ONE-SHOT, kept as the protocol of the rewrite (HANDOVER §11.8 L6). Guarded:
refuses to run when no chain is left. Run from the toolbox root.

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5 -- R4 Level L6, WP-64b.
"""
import re
import sys

FILES = ['tests/test_adjustVCE.m', 'tests/test_ellmer.m', 'tests/test_fit.m',
         'examples/demo_grace.m', 'examples/demo_messaging.m',
         'examples/demo_workflow.m']
CONF = {'groupComponentsBy': 'GroupBy', 'stochasticModel': 'Shape',
        'varianceComponents': 'Scale', 'iteration': 'Outer', 'joint': 'Joint',
        'componentsPerEpoch': 'Axes'}
CHAIN = re.compile(r'\.(' + '|'.join(CONF) + r')\(')


def split_statements(lines):
    """Join `...` continuations; yield (first, last, text, comment)."""
    out = []
    buf, start, comment = '', None, ''
    for i, ln in enumerate(lines):
        stripped = ln.strip()
        if buf == '':
            start = i
        code, _, cmt = ln.partition('%') if not stripped.startswith('%') else ('', '', ln)
        if '...' in code:
            buf += code.split('...')[0].rstrip() + ' '
            if cmt.strip():
                comment += ' ' + cmt.strip()
            continue
        buf += ln if not code else code.rstrip()
        if cmt.strip() and code:
            comment += ' ' + cmt.strip()
        out.append((start, i, buf, comment.strip()))
        buf, comment = '', ''
    return out


def match_paren(s, i):
    """Index of the ')' matching the '(' at s[i]."""
    depth = 0
    q = None
    k = i
    while k < len(s):
        c = s[k]
        if q:
            if c == q:
                q = None
        elif c in "'\"":
            # a quote starts a string unless it is a transpose (after word/paren)
            if c == "'" and k > 0 and (s[k-1].isalnum() or s[k-1] in ")]}'."):
                pass
            else:
                q = c
        elif c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return k
        k += 1
    raise ValueError('unbalanced: ' + s)


def parse_chain(text):
    """Split `lhs = base.m1(a1).m2(a2)...;` -> lhs, base, [(name, args)]."""
    m = re.match(r'^(\s*)(?:([\w\[\], ~]+?)\s*=\s*)?(.*?);\s*$', text)
    if not m:
        return None
    indent, lhs, rest = m.group(1), m.group(2), m.group(3)
    # base: either a variable or an inline call ending with ')'
    if rest.startswith('lsqvce.engine.adjustVCE('):
        j = match_paren(rest, rest.index('('))
        base, chain = rest[:j + 1], rest[j + 1:]
    else:
        mm = re.match(r'^(\w+)(.*)$', rest)
        base, chain = mm.group(1), mm.group(2)
    calls = []
    k = 0
    chain = chain.strip()
    while k < len(chain):
        while k < len(chain) and chain[k].isspace():
            k += 1
        if k >= len(chain):
            break
        mm = re.match(r'\.(\w+)\(', chain[k:])
        if not mm:
            raise ValueError('chain parse: ' + chain[k:])
        name = mm.group(1)
        po = k + len(mm.group(0)) - 1
        pc = match_paren(chain, po)
        calls.append((name, chain[po + 1:pc].strip()))
        k = pc + 1
    return indent, lhs, base, calls


def groups_to_nv(calls):
    nv, rest = [], []
    for name, args in calls:
        if name not in CONF:
            rest.append((name, args))
            continue
        key = CONF[name]
        if key == 'GroupBy':
            nv.append("'GroupBy', %s" % args)
        elif key == 'Axes':
            nv.append("'Axes', %s" % args)
        else:
            if key == 'Joint':
                nv.append("'VceMode', 'joint'")
            nv.append("'%s', struct(%s)" % (key, args) if args else None)
    return [x for x in nv if x], rest


def insert_nv(call_text, nv):
    """Insert NV pairs before the final ')' of a constructor/builder call."""
    j = call_text.rindex(')')
    return call_text[:j] + ', ' + ', '.join(nv) + call_text[j:]


def rewrite_file(path):
    src = open(path, encoding='utf-8').read()
    lines = src.split('\n')
    stmts = split_statements(lines)
    edits = {}     # first line index -> (last, new text)
    deleted = set()
    n = 0
    for idx, (a, b, text, cmt) in enumerate(stmts):
        if text.strip().startswith('%') or not CHAIN.search(text):
            continue
        parsed = parse_chain(text.strip('\n'))
        if parsed is None:
            raise ValueError('cannot parse %s:%d' % (path, a + 1))
        indent, lhs, base, calls = parsed
        nv, rest = groups_to_nv(calls)
        tail = ''.join('.%s(%s)' % (nm, ar) for nm, ar in rest)
        trailing = ('   % ' + cmt) if cmt else ''
        if base.startswith('lsqvce.engine.adjustVCE('):
            new = '%s%s%s%s;%s' % (indent, (lhs + ' = ') if lhs else '', insert_nv(base, nv), tail, trailing)
            edits[a] = (b, new)
        else:
            # find the builder statement of `base` before this one
            found = None
            for k in range(idx - 1, -1, -1):
                ta = stmts[k][2]
                if re.match(r'^\s*(\[\s*%s\b[^\]]*\]|%s)\s*=\s*.*(makeVCEBuilder|lsqvce\.engine\.adjustVCE)\(' % (base, base), ta):
                    found = k
                    break
                if re.match(r'^\s*function\b', ta):
                    break
            if found is None:
                raise ValueError('no builder for %s at %s:%d' % (base, path, a + 1))
            fa, fb, ftext, fcmt = stmts[found]
            ftext2 = insert_nv(ftext.rstrip(), nv)
            if fcmt:
                ftext2 += '   % ' + fcmt
            edits[fa] = (fb, ftext2)
            if tail == '' and (lhs is None or lhs.strip() == base):
                deleted.add(a)
                for r in range(a, b + 1):
                    deleted.add(r)
            else:
                new = '%s%s%s%s;%s' % (indent, (lhs + ' = ') if lhs else '', base, tail, trailing)
                edits[a] = (b, new)
        n += 1
    # apply
    out = []
    i = 0
    while i < len(lines):
        if i in edits:
            last, new = edits[i]
            out.append(new)
            i = last + 1
        elif i in deleted:
            i += 1
        else:
            out.append(lines[i])
            i += 1
    open(path, 'w', encoding='utf-8').write('\n'.join(out))
    return n


def main():
    total = 0
    for f in FILES:
        total += rewrite_file(f)
    if total == 0:
        print('nothing to rewrite (already done) -- refusing to run twice')
        sys.exit(1)
    print('rewrote %d chain statements' % total)


if __name__ == '__main__':
    main()
