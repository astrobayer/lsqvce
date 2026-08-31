#!/usr/bin/env python3
"""Route every hardcoded banner literal through lsqvce.util.verbosity.banner.

WHY
    The banner vocabulary was centralised in verbosity.banner so every log line
    aligns, but only adjustVCE and BlockAccumulator's loop were converted. The
    rest of the toolbox still wrote its own literals, and the MATLAB log proved
    it: '[resid ]' appeared 499 times and '[resid   ]' 63 times -- the SAME tag
    at two widths, in one run. This rewrites those literals to '%s' and injects
    banner('tag') as the first argument, so the width is decided in exactly one
    place.

WHAT IT HANDLES
    Two call shapes, distinguished by where the format string ends:
        log.msg('  [tag] text %d\n', a)             -> single quoted string
        log.msg(['  [tag] text ' ... 'more\n'], a)  -> bracketed concatenation
    In both cases the banner argument must land immediately AFTER the complete
    format string, which means finding the matching close of the bracket, not
    just the next quote.

USAGE
    python3 dev/docs_src/convert_banners.py            # from the toolbox root

Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""

import re
import sys

BANNER = 'lsqvce.util.verbosity.banner'

# (path, literal-as-written, tag)
TARGETS = [
    ('+lsqvce/residuals.m',          '[resid ]',   'resid'),
    ('+lsqvce/residuals.m',          '[resid   ]', 'resid'),
    ('+lsqvce/redundancy.m',         '[redund  ]', 'redund'),
    ('+lsqvce/vce.m',                '[vce     ]', 'vce'),
    ('+lsqvce/streamContribution.m', '[contrib]',  'contrib'),
    ('+lsqvce/poolGroups.m',         '[pool    ]', 'pool'),
    ('+lsqvce/selectArOrder.m',      '[selectAR]', 'selectAR'),
    ('+lsqvce/outliers.m',           '[outlier ]', 'outlier'),
    ('+lsqvce/mapreduce.m',          '[mapred  ]', 'mapred'),
    ('+lsqvce/iterateNoiseModel.m',  '[psdIt   ]', 'psdIt'),
    ('+lsqvce/contribution.m',       '[contrib ]', 'contrib'),
    ('+lsqvce/adjustVCEjoint.m',     '[joint]',    'joint'),
    ('+lsqvce/BlockAccumulator.m',   '[solve   ]', 'solve'),
    ('+lsqvce/BlockAccumulator.m',   '[addSys  ]', 'addSys'),
]


def end_of_format(text, start):
    """Index just past the format string beginning at `start` in `text`.

    `text` is the whole file with MATLAB line continuations already joined, so a
    bracketed format sits on one logical line. Returns None if the shape is not
    recognised, so the caller can skip rather than corrupt the source.
    """
    i = start
    while i < len(text) and text[i] in ' \t':
        i += 1
    if i >= len(text):
        return None

    if text[i] == '[':
        depth = 0
        instr = False
        while i < len(text):
            c = text[i]
            if instr:
                if c == "'":
                    if i + 1 < len(text) and text[i + 1] == "'":
                        i += 2
                        continue
                    instr = False
            elif c == "'":
                instr = True
            elif c == '[':
                depth += 1
            elif c == ']':
                depth -= 1
                if depth == 0:
                    return i + 1
            i += 1
        return None

    if text[i] == "'":
        i += 1
        while i < len(text):
            c = text[i]
            if c == "'":
                if i + 1 < len(text) and text[i + 1] == "'":
                    i += 2
                    continue
                return i + 1
            i += 1
        return None

    return None


def convert_file(path, edits):
    src = open(path, encoding='utf-8').read()
    # Join continuations so a bracketed format is one logical line, but REMEMBER
    # the original text: we rewrite the original, using the joined copy only to
    # locate the end of each format string.
    changed = 0
    for lit, tag in edits:
        while True:
            m = re.search(r'log\.msg\(', src)
            hit = None
            for m in re.finditer(r'log\.msg\(', src):
                seg_end = src.find('\n\n', m.end())
                seg = src[m.end(): seg_end if seg_end > 0 else len(src)]
                if lit in seg[:400]:
                    hit = m
                    break
            if hit is None:
                break

            # Join continuations from the call start so bracket matching works.
            tail = src[hit.end():]
            joined = re.sub(r"\.\.\.[^\n]*\n\s*", " ", tail)
            e = end_of_format(joined, 0)
            if e is None:
                print('  ! could not parse the format string in %s' % path)
                break

            # Map the joined offset back to an offset in the raw tail by walking
            # both, skipping the continuation runs in the raw text.
            ri = ji = 0
            while ji < e and ri < len(tail):
                mcont = re.match(r"\.\.\.[^\n]*\n\s*", tail[ri:])
                if mcont:
                    ri += mcont.end()
                    ji += 1          # the run collapsed to a single space
                    continue
                ri += 1
                ji += 1
            raw_end = hit.end() + ri

            before = src[:raw_end]
            after = src[raw_end:]
            if lit not in before:
                print('  ! literal moved unexpectedly in %s' % path)
                break
            before = before.replace(lit, '%s', 1)
            src = before + ", " + BANNER + "('" + tag + "')" + after
            changed += 1
    open(path, 'w', encoding='utf-8').write(src)
    return changed


def main():
    byfile = {}
    for path, lit, tag in TARGETS:
        byfile.setdefault(path, []).append((lit, tag))
    total = 0
    for path, edits in byfile.items():
        n = convert_file(path, edits)
        total += n
        print('  %-35s %d site(s)' % (path, n))
    print('\n%d site(s) converted.' % total)
    return 0


if __name__ == '__main__':
    sys.exit(main())
