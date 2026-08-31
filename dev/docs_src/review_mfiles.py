#!/usr/bin/env python3
"""Review every .m file for documentation completeness and robustness.

This is a REVIEW tool, not a gate: it reports what is missing so a human can
judge, and it deliberately does not fail the build. audit_tree.py is the gate.

Run from the toolbox root:  python3 dev/docs_src/review_mfiles.py

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
import glob
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_tree import code_lines            # the validated comment stripper

FINDINGS = defaultdict(list)


def note(kind, msg):
    FINDINGS[kind].append(msg)


def header_block(path):
    """The contiguous %-comment block immediately after the first function/classdef."""
    out, started = [], False
    for raw in open(path, encoding='utf-8', errors='replace'):
        s = raw.rstrip('\n')
        if not started:
            if re.match(r'\s*(function|classdef)\b', s):
                started = True
            continue
        if s.lstrip().startswith('%'):
            out.append(s.lstrip()[1:])
        elif s.strip() == '':
            if out:                      # a blank line ends the header block
                break
        else:
            break
    return '\n'.join(out)


def public_functions(path):
    """Top-level function/classdef names declared in a file."""
    names = []
    for _, s in code_lines(path):
        m = re.match(r'\s*(?:function\s+(?:\[[^\]]*\]|[\w~, ]+)\s*=\s*|function\s+|classdef\s+)(\w+)', s)
        if m:
            names.append(m.group(1))
    return names


def main():
    if not os.path.isdir('+lsqvce'):
        sys.exit('run from the toolbox root')

    files = sorted(glob.glob('+lsqvce/**/*.m', recursive=True))
    print('reviewing %d files in +lsqvce\n' % len(files))
    print('=' * 66)

    for f in files:
        base = os.path.basename(f)
        if base == 'Contents.m':
            continue
        h = header_block(f)
        name = os.path.splitext(base)[0]

        # --- 1. a header at all -------------------------------------------
        if not h.strip():
            note('no-help', f)
            continue

        # --- 2. H1 line names the function --------------------------------
        # A package function's H1 is conventionally LSQVCE.NAME or
        # LSQVCE.UTIL.NAME, not the bare name -- that is correct MATLAB style,
        # so compare against the LAST dotted component.
        first = h.strip().split('\n')[0].strip()
        h1word = first.split()[0] if first.split() else ''
        if h1word.split('.')[-1].upper() != name.upper():
            note('h1-mismatch', '%s: H1 starts %r' % (f, first[:44]))

        # --- 3. INPUT / OUTPUT sections -----------------------------------
        # A classdef documents properties/methods instead; a function that
        # takes or returns anything should say what.
        is_class = bool(re.search(r'^\s*classdef\b', open(f, encoding='utf-8',
                                                          errors='replace').read(), re.M))
        if not is_class:
            sig = ''
            for _, s in code_lines(f):
                if re.match(r'\s*function\b', s):
                    sig = s
                    break
            takes = '(' in sig and not re.search(r'\(\s*\)', sig)
            gives = '=' in sig.split('function', 1)[-1].split('(')[0]
            # Accept the variants actually in use: 'INPUT', 'INPUTS',
            # 'OUTPUT V', 'OUTPUT REJ', and 'OPTIONS' for a name-value-only
            # front door. Being stricter than the house style just reports
            # correct files as defective, which is worse than not checking.
            if takes and not re.search(r'^\s*(INPUTS?\b|OPTIONS\b|ARGUMENTS\b)',
                                       h, re.M | re.I):
                note('no-input-doc', f)
            if gives and not re.search(r'^\s*OUTPUTS?\b', h, re.M | re.I):
                note('no-output-doc', f)

        # --- 4. provenance ------------------------------------------------
        if not re.search(r'Claude|Anthropic', h):
            note('no-provenance', f)

        # --- 5. errors thrown but not documented --------------------------
        src = open(f, encoding='utf-8', errors='replace').read()
        ids = set(re.findall(r"error\(\s*'(lsqvce[:\w]+)'", src))
        if ids and not re.search(r'^\s*ERRORS?\b', h, re.M | re.I):
            note('undocumented-errors',
                 '%s: throws %s' % (f, ', '.join(sorted(ids)[:3])))

        # --- 6. See also --------------------------------------------------
        if not re.search(r'See also', h, re.I):
            note('no-see-also', f)

        # --- 7. long functions without section banners --------------------
        n = sum(1 for _ in code_lines(f))
        if n > 400 and src.count('% ---') < 3:
            note('long-unsectioned', '%s (%d code lines)' % (f, n))

    # --- report ------------------------------------------------------------
    order = ['no-help', 'h1-mismatch', 'no-input-doc', 'no-output-doc',
             'undocumented-errors', 'no-provenance', 'no-see-also',
             'long-unsectioned']
    titles = {
        'no-help': 'NO HELP BLOCK AT ALL',
        'h1-mismatch': 'H1 line does not start with the function name',
        'no-input-doc': 'takes arguments but has no INPUT section',
        'no-output-doc': 'returns values but has no OUTPUT section',
        'undocumented-errors': 'throws lsqvce:* errors with no ERRORS section',
        'no-provenance': 'no provenance line',
        'no-see-also': 'no "See also"',
        'long-unsectioned': 'long file with few section banners',
    }
    total = 0
    for k in order:
        v = FINDINGS.get(k, [])
        if not v:
            continue
        total += len(v)
        print('\n%s  (%d)' % (titles[k], len(v)))
        for item in v[:14]:
            print('   ', item)
        if len(v) > 14:
            print('    ... and %d more' % (len(v) - 14))

    print('\n' + '=' * 66)
    print('%d observation(s). This is a REVIEW, not a gate -- judge each one.'
          % total)


if __name__ == '__main__':
    main()
