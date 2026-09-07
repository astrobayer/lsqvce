#!/usr/bin/env python3
"""Fault-injection self-test of the R5 M0 audit checks (16, 19, 20, 21).

WHY THIS EXISTS
    "A check without a fixture proving it fires is not a check"
    (VALIDATION_GUIDE Part 3): twice in the source project a static audit's
    silence was read as evidence about something it never looked at. Every
    check added in M0 therefore ships with a fault it must catch. This script
    copies the tree to a scratch folder, injects ONE fault per case, runs the
    relevant check function of audit_tree.py in that folder and asserts that
    (a) the healthy copy is silent and (b) the injected copy fires with the
    expected finding kind. It never touches the real tree.

USAGE
    python3 dev/docs_src/audit_selftest.py        # from the toolbox root
    echo $?                                        # 0 = every fault caught

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5.1 -- R5 Level M0, WP-78.
"""
import importlib.util
import os
import shutil
import sys
import tempfile

ROOT = os.getcwd()
if not os.path.isdir('+lsqvce'):
    sys.exit('run from the toolbox root (the folder holding +lsqvce)')

spec = importlib.util.spec_from_file_location(
    'audit_tree', os.path.join(ROOT, 'dev', 'docs_src', 'audit_tree.py'))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

COPY = ['+lsqvce', 'html/index.html', 'helptoc.xml', 'dev/docs_src/tiers.csv', 'dev/docs_src/catch_baseline.csv',
        'dev/docs_src/error_ids.csv', 'dev/docs_src/layering_baseline.csv',
        'dev/tests/category_exemptions.csv', 'dev/tests']


def scratch():
    d = tempfile.mkdtemp(prefix='lsqvce_selftest_')
    for item in COPY:
        src = os.path.join(ROOT, item)
        dst = os.path.join(d, item)
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
    return d


def run(check):
    audit.FINDINGS.clear()
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        check()
    return list(audit.FINDINGS)


def append(path, text):
    with open(path, 'a', encoding='utf-8') as fh:
        fh.write(text)


CASES = []


def case(name, check, kind, mutate, expect_fire=True):
    CASES.append((name, check, kind, mutate, expect_fire))


# ---- 16: seven categories -------------------------------------------------
case('16 metamorphic is legal vocabulary', audit.check_test_categories, 'cat',
     lambda d: append(os.path.join(d, 'dev/tests/test_zzz_selftest.m'),
                      'function tests = test_zzz_selftest\n'
                      '    tests = functiontests(localfunctions);\nend\n'
                      'function testMeta(tc)\n% categories: metamorphic | tags: concurrency\n'
                      '    lsqvce.fit;\nend\n'),
     expect_fire=False)
case('16 unknown category fires', audit.check_test_categories, 'cat',
     lambda d: append(os.path.join(d, 'dev/tests/test_zzz_selftest.m'),
                      'function tests = test_zzz_selftest\n'
                      '    tests = functiontests(localfunctions);\nend\n'
                      'function testBogus(tc)\n% categories: bogus\n'
                      '    lsqvce.fit;\nend\n'))

# ---- 19: catch ratchet ----------------------------------------------------
case('19 new bare catch fires', audit.check_catch_ratchet, 'catch',
     lambda d: append(os.path.join(d, '+lsqvce/+util/getfielddef.m'),
                      '\nfunction zzzSelftest()\n    try\n        x = 1;\n'
                      '    catch\n        x = 2;\n    end\nend\n'))
case('19 stale baseline row fires', audit.check_catch_ratchet, 'catch',
     lambda d: append(os.path.join(d, 'dev/docs_src/catch_baseline.csv'),
                      'engine.zzzNoSuchUnit:1,bare\n'))

# ---- 20: identifier catalogue --------------------------------------------
case('20 uncatalogued id fires', audit.check_id_catalogue, 'id-catalogue',
     lambda d: append(os.path.join(d, '+lsqvce/+util/getfielddef.m'),
                      "\nfunction zzzSelftest()\n"
                      "    error('lsqvce:selftest:injected', 'x');\nend\n"))
case('20 stale catalogue row fires', audit.check_id_catalogue, 'id-catalogue',
     lambda d: append(os.path.join(d, 'dev/docs_src/error_ids.csv'),
                      'lsqvce:selftest:ghost,util.getfielddef,error,v4.0\n'))

# ---- 21: layering ratchet -------------------------------------------------
case('21 new upward call fires', audit.check_layering_ratchet, 'layering',
     lambda d: append(os.path.join(d, '+lsqvce/+kernel/dctI.m'),
                      '\nfunction zzzSelftest()\n    lsqvce.fit(struct());\nend\n'))
case('21 stale baseline row fires', audit.check_layering_ratchet, 'layering',
     lambda d: append(os.path.join(d, 'dev/docs_src/layering_baseline.csv'),
                      'kernel.dctI,report,KERNEL,USER\n'))



# ---- E-12 (R6 WP-93): the audit reads Windows paths like POSIX ones -------
# glob returns backslash-separated paths on Windows; every unit name, tier and
# index key is derived by '/'-string work. Simulate Windows on any platform by
# wrapping glob to hand back backslashes: the healthy tree must stay SILENT on
# checks 16 and 21 (before pglob: ~590 false category gaps and every layering
# baseline row reported stale).
def on_windows_paths(check):
    def wrapped():
        real = audit.glob.glob
        audit.glob.glob = lambda pat, **kw: [p.replace('/', '\\') for p in real(pat, **kw)]
        try:
            check()
        finally:
            audit.glob.glob = real
    wrapped.__name__ = check.__name__ + '_on_windows_paths'
    return wrapped


case('16 silent on Windows-style paths (E-12)', on_windows_paths(audit.check_test_categories),
     'cat', lambda d: None, expect_fire=False)
case('21 silent on Windows-style paths (E-12)', on_windows_paths(audit.check_layering_ratchet),
     'layering', lambda d: None, expect_fire=False)
case('21 new upward call still fires on Windows-style paths', on_windows_paths(audit.check_layering_ratchet),
     'layering',
     lambda d: append(os.path.join(d, '+lsqvce/+kernel/dctI.m'),
                      '\nfunction zzzSelftest()\n    lsqvce.fit(struct());\nend\n'))


# ---- 24 / 25 (R6 WP-94): help examples and the index ----------------------
case('24 example calling a deleted unit fires', audit.check_help_examples_call_live_units,
     'example',
     lambda d: append(os.path.join(d, '+lsqvce/+util/getfielddef.m'),
                      '\n%   EXAMPLE\n%     R = lsqvce.engine.zzzGone(1);\n'))
case('24 example calling a live unit is silent', audit.check_help_examples_call_live_units,
     'example',
     lambda d: append(os.path.join(d, '+lsqvce/+util/getfielddef.m'),
                      '\n%   EXAMPLE\n%     R = lsqvce.fit(data);\n'),
     expect_fire=False)
case('25 unit missing from index and helptoc fires', audit.check_index_lists_every_unit,
     'index',
     lambda d: append(os.path.join(d, '+lsqvce/+util/zzzNewUnit.m'),
                      'function zzzNewUnit()\n%ZZZNEWUNIT  x\nend\n'))
def main():
    failures = 0
    print('audit self-test: %d cases' % len(CASES))
    for name, check, kind, mutate, expect_fire in CASES:
        d = scratch()
        try:
            os.chdir(d)
            healthy = [f for f in run(check) if f[0] == kind]
            mutate(d)
            fired = [f for f in run(check) if f[0] == kind]
        finally:
            os.chdir(ROOT)
            shutil.rmtree(d, ignore_errors=True)
        ok_healthy = not healthy
        ok_fault = bool(fired) == expect_fire
        status = 'ok' if (ok_healthy and ok_fault) else 'FAIL'
        if status == 'FAIL':
            failures += 1
        print('  [%s] %-40s healthy:%d  after fault:%d  (expected fire=%s)'
              % (status, name, len(healthy), len(fired), expect_fire))
    if failures:
        print('%d self-test case(s) FAILED' % failures)
        return 1
    print('every fault caught, healthy copies silent')
    return 0


if __name__ == '__main__':
    sys.exit(main())
