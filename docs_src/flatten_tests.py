#!/usr/bin/env python3
"""Flatten tests/accum into tests/ and rename the release-named test files.

WHY (HANDOVER 1.7b)
    `tests/` held 7 files and `tests/accum/` held 12 with no rule saying which
    went where. Three problems:
      1. THE SPLIT IS A TRAP. runtests('tests') does not recurse, so the obvious
         command silently ran 7 files and skipped 12. The handover warned about
         it twice; a warning is a workaround, not a fix.
      2. 'accum' stopped describing its contents (it held the VCE loop tests,
         contribution, precision, reliability -- none about accumulation).
      3. Five files were named after RELEASES, not subjects: test_v3/v4/v7/v8/v9
         meant "everything added in version N". Their headers show they ARE
         coherent, just mislabelled. test_recover was named after a function
         REMOVED on 2026-07-18 (recover -> adjustBlock).

WHAT IT DOES
    Moves every tests/accum/*.m up one level, renames the six misnamed files,
    rewrites the `function tests = <name>` line inside each (MATLAB requires it
    to match the file name), fixes the one path computation that assumed the
    extra directory level, and updates every reference in the runners, the
    setup script, README.md, HANDOVER.md and html/index.html.

    Idempotent-ish: it reports and skips anything already moved.

USAGE
    python3 docs_src/flatten_tests.py        # from the toolbox root

Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""

import os
import re
import shutil
import sys

# old stem -> new stem. Only the misnamed ones are renamed; thematically named
# files keep their names so the churn stays as small as the problem allows.
RENAME = {
    'test_v3':      'test_vceMath',       # leverage, Helmert VCE, snooping, reg, multi-axis
    'test_v4':      'test_parallel',      # mapreduce, addSystems, redundancy, PSD iteration
    'test_v7':      'test_vceFunction',   # lsqvce.engine.vce contract + diagnoseAccum bridge
    'test_v8':      'test_signalCore',    # outliers, varWhiten, ar2arma
    'test_v9':      'test_streamUtils',   # selectArAcv, whitenessCorr, collectResiduals, ...
    'test_recover': 'test_adjustBlock',   # 'recover' was removed 2026-07-18
}

TEXT_FILES = ['tests/run_all_lsqvce.m', 'tests/run_lsqvce.m', 'lsqvce_setup.m',
              'README.md', 'HANDOVER.md', 'html/index.html']


def move_and_rename():
    src_dir = os.path.join('tests', 'accum')
    if not os.path.isdir(src_dir):
        print('  tests/accum is already gone; skipping the move')
        return []
    moved = []
    for fn in sorted(os.listdir(src_dir)):
        if not fn.endswith('.m'):
            continue
        stem = fn[:-2]
        new_stem = RENAME.get(stem, stem)
        src = os.path.join(src_dir, fn)
        dst = os.path.join('tests', new_stem + '.m')
        if os.path.exists(dst):
            print('  ! %s already exists, skipping %s' % (dst, src))
            continue
        shutil.move(src, dst)
        moved.append((stem, new_stem, dst))
        print('  %-18s -> tests/%s.m' % ('accum/' + fn, new_stem))
    try:
        os.rmdir(src_dir)
        print('  removed empty tests/accum/')
    except OSError as e:
        print('  ! could not remove tests/accum: %s' % e)
    return moved


def fix_entry_points(moved):
    """MATLAB requires the top function name to match the file name."""
    for old_stem, new_stem, path in moved:
        if old_stem == new_stem:
            continue
        s = open(path, encoding='utf-8').read()
        pat = re.compile(r'^(function\s+tests\s*=\s*)%s\b' % re.escape(old_stem),
                         re.M)
        s2, n = pat.subn(r'\g<1>' + new_stem, s)
        if n == 0:
            print('  ! no entry point "function tests = %s" in %s'
                  % (old_stem, path))
        # The H1 line and any self-reference in the help block.
        s2 = re.sub(r'\b%s\b' % re.escape(old_stem.upper()),
                    new_stem.upper(), s2)
        open(path, 'w', encoding='utf-8').write(s2)
        print('  entry point + H1 updated in %s' % path)


def fix_paths():
    """One test walked up two levels to reach examples/; now it is one."""
    p = os.path.join('tests', 'test_separable.m')
    if not os.path.exists(p):
        return
    s = open(p, encoding='utf-8').read()
    old = ("    here = fileparts(mfilename('fullpath'));               % tests/accum\n"
           "    exDir = fullfile(here, '..', '..', 'examples');")
    new = ("    here = fileparts(mfilename('fullpath'));               % tests/\n"
           "    exDir = fullfile(here, '..', 'examples');")
    if old in s:
        open(p, 'w', encoding='utf-8').write(s.replace(old, new, 1))
        print('  fixed the examples/ path in tests/test_separable.m')
    else:
        print('  ! path computation in test_separable.m not in the expected form')


def fix_references(moved):
    ren = {o: n for o, n, _ in moved if o != n}
    for f in TEXT_FILES:
        if not os.path.exists(f):
            continue
        s = open(f, encoding='utf-8').read()
        o = s
        # Path forms first, so the stem rename below cannot corrupt them.
        s = s.replace("fullfile(root,'tests','accum',", "fullfile(root,'tests',")
        s = s.replace("fullfile(root, 'tests', 'accum',", "fullfile(root, 'tests',")
        s = s.replace("fullfile(root,'tests','accum')", "fullfile(root,'tests')")
        s = s.replace('tests/accum/', 'tests/')
        s = s.replace('tests\\accum\\', 'tests\\')
        s = s.replace('tests/accum', 'tests')
        for old_stem, new_stem in ren.items():
            s = re.sub(r'\b%s\b' % re.escape(old_stem), new_stem, s)
        if s != o:
            open(f, 'w', encoding='utf-8').write(s)
            print('  updated references in %s' % f)

    # Cross-references inside the test files themselves.
    for fn in sorted(os.listdir('tests')):
        if not fn.endswith('.m'):
            continue
        p = os.path.join('tests', fn)
        s = open(p, encoding='utf-8').read()
        o = s
        s = s.replace('tests/accum/', 'tests/').replace('tests/accum', 'tests')
        for old_stem, new_stem in ren.items():
            s = re.sub(r'\b%s\b' % re.escape(old_stem), new_stem, s)
        if s != o:
            open(p, 'w', encoding='utf-8').write(s)
            print('  updated cross-references in %s' % p)


def main():
    if not os.path.isdir('tests'):
        sys.exit('run from the toolbox root (the folder holding tests/)')
    print('flattening tests/')
    moved = move_and_rename()
    fix_entry_points(moved)
    fix_paths()
    fix_references(moved)
    print('\ndone. Run docs_src/audit_tree.py to confirm registration.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
