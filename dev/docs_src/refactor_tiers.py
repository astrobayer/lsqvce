#!/usr/bin/env python3
"""R4-D1 (WP-64a): move the ENGINE and KERNEL tiers into +lsqvce/+engine and
+lsqvce/+kernel, by the ONE tier authority dev/docs_src/tiers.csv.

ONE-SHOT, MECHANICAL, RUN ONCE from the toolbox root on 2026-09-01. Kept in
the tree as the protocol of what was done; after the run it is guarded and
refuses to run again (the tiers.csv it needs no longer exists in the input
shape). No hand edits are part of this step -- every change below is scripted
so it can be audited and, if needed, reasoned about as a whole.

What it does
  1. git mv every ENGINE row to +lsqvce/+engine/<name>.m, every KERNEL row
     (not util.*) to +lsqvce/+kernel/<name>.m; USER rows and +util stay.
  2. Private helpers: the six engine-only ones to +lsqvce/+engine/private/;
     runMeta -- used by USER (Results, fit) AND ENGINE (adjustBlock,
     adjustVCE) -- becomes the public engine helper lsqvce.engine.runMeta.
     (r1 correction: probeStream IS called, by the private probeMatrix; it
     was restored to +engine/private. The promote regex also rewrote runMeta's
     own declaration line -- fixed by hand in r1; the guard below keeps the
     script from ever running again.)
  3. Rewrites every qualified reference lsqvce.<name> / LSQVCE.<NAME> to its
     new home in code, tests, examples, tooling, HTML, helptoc, README,
     DATA_CONTRACT, SETUP_GUIDE, CLAUDE.md. Error/warning IDENTIFIERS
     (lsqvce:<unit>:<rule>) are API and are NOT touched. HANDOVER.md and
     HISTORY_v1.md are history and are NOT rewritten.
  4. Renames html/lsqvce.<name>.html to lsqvce.<tier>.<name>.html.
  5. Regenerates +lsqvce/Contents.m (USER tier + package pointers) and writes
     +lsqvce/+engine/Contents.m, +lsqvce/+kernel/Contents.m from tiers.csv.
  6. Rewrites tiers.csv names to <tier>.<name> and the unit column of
     tests/category_exemptions.csv.
  7. Stamps every moved file with one help line naming the move.

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5 -- R4 Level L5, WP-64a.
"""
import csv
import os
import re
import subprocess
import sys

_STALE_GUARD = True   # set False ONLY to re-run on a tree in the pre-move shape

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # dev/docs_src -> root (L8)
os.chdir(ROOT)

TIER_DIR = {'ENGINE': '+engine', 'KERNEL': '+kernel'}
TIER_PKG = {'ENGINE': 'engine', 'KERNEL': 'kernel'}
ENGINE_PRIVATE = ['applyCovOverride', 'arcBatchWhiten', 'jointArcParts',
                  'makeKeyFcn', 'probeMatrix']
PROMOTE = ['runMeta']
DELETE = []   # probeStream is called by the private probeMatrix (missed by a
              # root-only grep on 2026-09-01); restored to +engine/private in r1
SKIP_FILES = {'HANDOVER.md', 'HISTORY_v1.md'}
TEXT_EXT = {'.m', '.py', '.html', '.xml', '.md', '.csv', '.txt'}


def sh(*args):
    subprocess.run(args, check=True)


def rows():
    return list(csv.DictReader(open('dev/docs_src/tiers.csv', encoding='utf-8')))


def main():
    R = rows()
    if _STALE_GUARD and any(r['name'].startswith(('engine.', 'kernel.')) for r in R):
        sys.exit('refactor_tiers.py: tree is already in the post-move shape; '
                 'this script ran once (2026-09-01) and is kept as protocol only.')
    moves = {}                       # name -> tier
    for r in R:
        if r['tier'] in TIER_DIR and not r['name'].startswith('util.'):
            moves[r['name']] = r['tier']
    for t in TIER_DIR.values():
        os.makedirs(os.path.join('+lsqvce', t), exist_ok=True)
    os.makedirs('+lsqvce/+engine/private', exist_ok=True)

    # 1. move the tier files
    for name, tier in sorted(moves.items()):
        sh('git', 'mv', '+lsqvce/%s.m' % name, '+lsqvce/%s/%s.m' % (TIER_DIR[tier], name))
        page = 'html/lsqvce.%s.html' % name
        if os.path.exists(page):
            sh('git', 'mv', page, 'html/lsqvce.%s.%s.html' % (TIER_PKG[tier], name))
    # 2. private helpers
    for n in ENGINE_PRIVATE:
        sh('git', 'mv', '+lsqvce/private/%s.m' % n, '+lsqvce/+engine/private/%s.m' % n)
    for n in PROMOTE:
        sh('git', 'mv', '+lsqvce/private/%s.m' % n, '+lsqvce/+engine/%s.m' % n)
    for n in DELETE:
        sh('git', 'rm', '-q', '+lsqvce/private/%s.m' % n)
    if os.path.isdir('+lsqvce/private') and not os.listdir('+lsqvce/private'):
        os.rmdir('+lsqvce/private')

    # 3. rewrite references
    rename = {n: 'lsqvce.%s.%s' % (TIER_PKG[t], n) for n, t in moves.items()}
    for n in PROMOTE:
        rename[n] = 'lsqvce.engine.%s' % n
    pat_lc = re.compile(r'\blsqvce\.(%s)\b' % '|'.join(map(re.escape, sorted(rename, key=len, reverse=True))))
    up = {n.upper(): rename[n].upper() for n in rename}
    pat_uc = re.compile(r'\bLSQVCE\.(%s)\b' % '|'.join(map(re.escape, sorted(up, key=len, reverse=True))))
    # unqualified private-call sites of the promoted helper (root + engine).
    # r1 lesson: this must NOT match the helper's own 'function ... = runMeta('
    # declaration -- exclude lines starting with 'function'.
    pat_promote = re.compile(r'^(?!\s*function\b)(.*?)(?<![\w.])(%s)\s*\(' % '|'.join(PROMOTE), re.M)

    touched = 0
    for dirpath, dirnames, files in os.walk('.'):
        dirnames[:] = [d for d in dirnames if d not in ('.git', '__pycache__', 'node_modules')]
        for fn in files:
            path = os.path.join(dirpath, fn)
            if fn in SKIP_FILES or os.path.splitext(fn)[1] not in TEXT_EXT:
                continue
            if fn == 'refactor_tiers.py':
                continue
            try:
                s = open(path, encoding='utf-8').read()
            except UnicodeDecodeError:
                continue
            s2 = pat_lc.sub(lambda m: rename[m.group(1)], s)
            s2 = pat_uc.sub(lambda m: up[m.group(1)], s2)
            if path.endswith('.m') and ('+lsqvce' in path):
                s2 = pat_promote.sub(lambda m: 'lsqvce.engine.%s(' % m.group(1), s2)
            if s2 != s:
                open(path, 'w', encoding='utf-8').write(s2)
                touched += 1
    print('rewrote references in %d files' % touched)

    # 7. stamp moved files
    stamp = ('%%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/%s (R4 decision R4-D1,\n'
             '%%   WP-64a): the %s tier leaves the USER surface; every qualified\n'
             '%%   reference was rewritten mechanically by dev/docs_src/refactor_tiers.py.\n')
    for name, tier in moves.items():
        stamp_file('+lsqvce/%s/%s.m' % (TIER_DIR[tier], name), stamp % (TIER_DIR[tier], tier))
    for n in PROMOTE:
        stamp_file('+lsqvce/+engine/%s.m' % n,
                   '%   Moved 2026-09-01 (CEST, UTC+02:00) from +lsqvce/private to the public\n'
                   '%   ENGINE tier (R4-D1, WP-64a): used by USER (Results, fit) and ENGINE\n'
                   '%   (adjustBlock, adjustVCE), so it cannot stay package-private.\n')

    # 5. Contents.m
    write_contents(R, moves)

    # 6. tiers.csv and exemptions
    with open('dev/docs_src/tiers.csv', 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'tier', 'role'])
        for r in R:
            n = r['name']
            if n in moves:
                n = '%s.%s' % (TIER_PKG[moves[n]], n)
            w.writerow([n, r['tier'], r['role']])
        w.writerow(['engine.runMeta', 'ENGINE', 'immutable provenance record for Results.meta (promoted from private, WP-64a)'])
    p = 'tests/category_exemptions.csv'
    if os.path.exists(p):
        out = []
        for ln in open(p, encoding='utf-8'):
            parts = ln.rstrip('\n').split(',', 2)
            if len(parts) == 3 and parts[0] in moves:
                parts[0] = '%s.%s' % (TIER_PKG[moves[parts[0]]], parts[0])
            out.append(','.join(parts) + '\n')
        open(p, 'w', encoding='utf-8').write(''.join(out))
    print('done: %d files moved, %d private relocated, %d promoted, %d deleted'
          % (len(moves), len(ENGINE_PRIVATE), len(PROMOTE), len(DELETE)))


def stamp_file(path, stamp):
    lines = open(path, encoding='utf-8').read().split('\n')
    i = 0
    while i < len(lines) and not re.match(r'^\s*(function|classdef)\b', lines[i]):
        i += 1
    j = i + 1
    while j < len(lines) and lines[j].lstrip().startswith('%'):
        j += 1
    lines[j:j] = stamp.rstrip('\n').split('\n')
    open(path, 'w', encoding='utf-8').write('\n'.join(lines))


def write_contents(R, moves):
    root_lines = open('+lsqvce/Contents.m', encoding='utf-8').read().split('\n')
    version_line = root_lines[0]
    user = [r for r in R if r['tier'] == 'USER']
    root = [version_line,
            '%',
            '% lsqvce -- least-squares adjustment with variance component estimation.',
            '%',
            '% THE USER TIER (R4, four names; HANDOVER §11.4). Start with the Getting',
            '% started page (doc lsqvce) and DATA_CONTRACT.md at the toolbox root.',
            '%']
    for r in user:
        root.append('%%   %-14s - %s' % (r['name'], r['role']))
    root += ['%',
             '% The engines and kernels below fit live in two sub-packages (R4-D1,',
             '% 2026-09-01); they are reachable for experts and never needed by a user:',
             '%   engine         - ENGINE tier: drivers, accumulators, loaders, VCE',
             '%                    machinery, diagnostics    -> help lsqvce.engine',
             '%   kernel         - KERNEL tier: pure numerical kernels (whitening,',
             '%                    spectra, transforms)      -> help lsqvce.kernel',
             '%   util           - shared utilities            -> help lsqvce.util',
             '%',
             '% Regenerated 2026-09-01 by dev/docs_src/refactor_tiers.py (WP-64a).',
             '']
    open('+lsqvce/Contents.m', 'w', encoding='utf-8').write('\n'.join(root))
    for tier, sub in (('ENGINE', '+engine'), ('KERNEL', '+kernel')):
        names = [r for r in R if r['tier'] == tier and not r['name'].startswith('util.')]
        c = ['%% lsqvce.%s -- the %s tier of lsqvce (R4-D1, WP-64a).' % (TIER_PKG[tier], tier),
             '%',
             '% Users reach these through lsqvce.fit / lsqvce.Options / lsqvce.Results;',
             '% experts call them directly. Roles from dev/docs_src/tiers.csv, the ONE tier',
             '% authority (audit check 18).',
             '%']
        for r in sorted(names, key=lambda r: r['name'].lower()):
            c.append('%%   %-20s - %s' % (r['name'], r['role']))
        if tier == 'ENGINE':
            c.append('%%   %-20s - %s' % ('runMeta', 'immutable provenance record for Results.meta'))
        c.append('')
        open('+lsqvce/%s/Contents.m' % sub, 'w', encoding='utf-8').write('\n'.join(c))


if __name__ == '__main__':
    main()
