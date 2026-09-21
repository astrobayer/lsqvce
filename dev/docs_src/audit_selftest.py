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
import re
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
        'dev/docs_src/error_ids.csv', 'dev/docs_src/layering_baseline.csv', 'dev/docs_src/superseded_claims.csv',
        'dev/tests/category_exemptions.csv', 'dev/tests', 'lsqvce_manual.pdf']


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

# ---- 27 (R6 WP-101): a documented option with no consumer ----------------
# MinPasses and LogFile shipped dead for a whole review; this is the ratchet.
# The scratch copy carries +lsqvce, so the fault appends a property to the
# REAL Options.m copy and the check must see it as unread.
def add_dead_prop(d):
    p = os.path.join(d, '+lsqvce/Options.m')
    s = open(p, encoding='utf-8').read()
    s = s.replace("        OnPass     = []",
                  "        OnPass     = []\n        ZzzSelftest (1,1) double = 1", 1)
    open(p, 'w', encoding='utf-8').write(s)


def add_used_prop(d):
    p = os.path.join(d, '+lsqvce/Options.m')
    s = open(p, encoding='utf-8').read()
    s = s.replace("        OnPass     = []",
                  "        OnPass     = []\n        ZzzSelftest (1,1) double = 1", 1)
    open(p, 'w', encoding='utf-8').write(s)
    q = os.path.join(d, '+lsqvce/fit.m')
    append(q, "\n% selftest: opts.ZzzSelftest is read here\n")


case('27 a dead Options property fires', audit.check_options_are_consumed, 'options',
     add_dead_prop)

# ---- 28 (R6 WP-106): superseded claims ------------------------------------
def add_help_claim(d):
    # INSIDE the help block (before See also), where the check reads claims
    p = os.path.join(d, '+lsqvce/fit.m')
    lines = open(p, encoding='utf-8').read().split('\n')
    lines.insert(2, '%   Robust: Not combinable with `Estimate="noise"` in one run')
    open(p, 'w', encoding='utf-8').write('\n'.join(lines))


case('28 a retired claim in a help block fires', audit.check_superseded_claims, 'superseded',
     add_help_claim)
case('28 the same phrase after See also (a stamp) is silent', audit.check_superseded_claims,
     'superseded',
     lambda d: append(os.path.join(d, '+lsqvce/fit.m'),
                      '\n%   See also NOTHING\n%   Updated: Not combinable with `Estimate="noise"` was retired\n'),
     expect_fire=False)
def add_manual_claim(d):
    os.makedirs(os.path.join(d, 'dev/docs_src/manual'), exist_ok=True)
    append(os.path.join(d, 'dev/docs_src/manual/03_data.md'),
           "\n```matlab\nD.arcIdx = {1:m};\n```\n")


case('28 a retired usage form in the Manual fires', audit.check_superseded_claims, 'superseded',
     add_manual_claim)
case('27 a property something reads is silent', audit.check_options_are_consumed, 'options',
     add_used_prop, expect_fire=False)

# ---- 26 (R6 WP-99): markdown pointers ------------------------------------
# The scratch copy carries no scope document at all, so the healthy run checks
# zero pointers and is silent by construction; each fault WRITES the document
# it needs. The two silent cases are not decoration: both forms cost a
# refinement when the check was measured against the real tree, where a
# backticked COMMAND ("python3 dev/docs_src/audit_tree.py") and a path written
# relative to its own document (dev/SETUP_GUIDE.md saying `tests/ci_gate.m`)
# were the only false positives it produced.
case('26 dead markdown link fires', audit.check_md_pointers, 'md-pointer',
     lambda d: append(os.path.join(d, 'README.md'),
                      '\nSee [the page](html/zzzNoSuchPage.html) for more.\n'))
case('26 dead backticked path fires', audit.check_md_pointers, 'md-pointer',
     lambda d: append(os.path.join(d, 'README.md'),
                      '\nRun `dev/docs_src/zzzNoSuchScript.py` first.\n'))
case('26 a command with a live path is silent', audit.check_md_pointers, 'md-pointer',
     lambda d: append(os.path.join(d, 'README.md'),
                      '\nRun `matlab -batch dev/tests/ci_gate.m` to gate.\n'),
     expect_fire=False)
case('26 a path relative to its document is silent', audit.check_md_pointers, 'md-pointer',
     lambda d: append(os.path.join(d, 'dev/DEVELOPMENT.md'),
                      '\nThe gate is `tests/ci_gate.m`.\n'),
     expect_fire=False)



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
# ---- 29 (R6 WP-109): the PDF Manual's version surface --------------------
# The defect this check was written for is the REAL one: the tree went to
# 6.1.1 in WP-108 and lsqvce_manual.pdf stayed at 6.1.0 because nothing read
# it. The fault here reproduces exactly that -- move the ONE version source
# ahead and leave the PDF where it is. The other two cases are the
# fails-closed half: an unreadable PDF and a missing one must be findings,
# because a check that goes quiet when its instrument breaks is worse than no
# check at all.
def bump_version(d):
    p = os.path.join(d, '+lsqvce/Contents.m')
    s = open(p, encoding='utf-8').read()
    open(p, 'w', encoding='utf-8').write(
        re.sub(r'^%(\s*)Version(\s+)\S+', r'%\g<1>Version\g<2>9.9.9', s, count=1, flags=re.M))


case('29 a stale PDF title fires', audit.check_manual_pdf_version, 'manual-pdf',
     bump_version)
case('29 an unreadable PDF fires', audit.check_manual_pdf_version, 'manual-pdf',
     lambda d: open(os.path.join(d, 'lsqvce_manual.pdf'), 'wb').write(b'%PDF-1.5\n%not a pdf\n'))
case('29 a missing PDF fires', audit.check_manual_pdf_version, 'manual-pdf',
     lambda d: os.remove(os.path.join(d, 'lsqvce_manual.pdf')))


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
