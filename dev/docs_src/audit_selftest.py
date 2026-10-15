#!/usr/bin/env python3
"""Fault-injection self-test of the audit checks (every check has a case since R7-32).

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

# R7-32: the scratch tree carries everything a check reads (the whole html
# set, the doc sources and builders, README / DEVELOPMENT, template), so
# every check -- not only the ten that had a case -- can be fault-injected.
COPY = ['+lsqvce', 'html', 'helptoc.xml', 'dev/docs_src', 'dev/tests', 'template',
        'examples', 'README.md', 'RELEASE_NOTES.md', 'dev/DEVELOPMENT.md', 'dev/HANDOVER.md',
        'dev/GRAVFIRE_INTERFACE.md',
        'lsqvce_manual.pdf', 'lsqvce_api.pdf']


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

# ---- 2, 3, 7, 9, 10, 15 (R7-26): the checks see dev/tests ----------------
# Since L8 these checks swept a top-level tests/ folder that no longer
# existed; a fault in a TEST file left the whole audit CLEAN. Each case
# injects into dev/tests, where the check had been blind.
def edit(path, old, new):
    s = open(path, encoding='utf-8').read()
    assert old in s, 'selftest anchor not found in %s' % path
    open(path, 'w', encoding='utf-8').write(s.replace(old, new, 1))


TFILE = 'dev/tests/test_util.m'
case('10 duplicate local function in a TEST file fires', audit.check_duplicate_functions,
     'duplicate', lambda d: append(os.path.join(d, TFILE),
                                   '\nfunction zzzDup()\nend\nfunction zzzDup()\nend\n'))
case('7 unbalanced end in a TEST file fires', audit.check_balance, 'balance',
     lambda d: append(os.path.join(d, TFILE), '\nfunction zzzOpen()\n    if true\nend\n'))
case("7 ['^' f(x) '$'] in a TEST file is silent", audit.check_balance, 'balance',
     lambda d: append(os.path.join(d, TFILE),
                      "\nfunction zzzCat()\n    p = ['^' upper('x') '$'];\nend\n"),
     expect_fire=False)
case('2 removed name in TEST code fires', audit.check_removed_names, 'names',
     lambda d: append(os.path.join(d, TFILE),
                      '\nfunction zzzOld()\n    covModelFactory(1);\nend\n'))
case('15 wrong H1 in a TEST file fires', audit.check_h1_names, 'h1',
     lambda d: edit(os.path.join(d, TFILE), '%TEST_UTIL ', '%TEST_ZZZ '))
case('3 nc as leading arg in TEST code fires', audit.check_signatures, 'signature',
     lambda d: append(os.path.join(d, TFILE),
                      '\nfunction zzzSig()\n    adjustBlock(nc, {}, 1);\nend\n'))
case('3 a missing adjustBlock fires', audit.check_signatures, 'signature',
     lambda d: os.remove(os.path.join(d, '+lsqvce/+engine/adjustBlock.m')))
case('9 a diverging block-spec literal fires', audit.check_spec_fields, 'spec',
     lambda d: edit(os.path.join(d, '+lsqvce/fit.m'),
                    "'dayList', dl, 'name', nm);", "'dayList', dl, 'name', nm, 'zzz', 1);"))
case('9 no block-spec literal at all fires', audit.check_spec_fields, 'spec',
     lambda d: [os.remove(os.path.join(d, p)) for p in
                ('+lsqvce/fit.m', '+lsqvce/+engine/private/normalizeSystems.m')])


# ---- 16 (R7-27): the matrix credits code in gated suites only ------------
def credit_by_comment_only(d):
    # strip every code reference to lsqvce.engine.resolveNc from the suites,
    # leaving a comment that names it: the contract cell must open
    for f in audit.pglob(os.path.join(d, 'dev/tests/test_*.m')):
        s = open(f, encoding='utf-8').read()
        if 'lsqvce.engine.resolveNc(' in s:
            s = s.replace('lsqvce.engine.resolveNc(', 'zzzGone(')
            open(f, 'w', encoding='utf-8').write(s)
    append(os.path.join(d, 'dev/tests/test_zzz_selftest.m'),
           'function tests = test_zzz_selftest\n'
           '    tests = functiontests(localfunctions);\nend\n'
           'function testOnlyAComment(tc)\n% categories: contract\n'
           '    % lsqvce.engine.resolveNc is named here, never called\n'
           '    tc.verifyTrue(true);\nend\n')


def credit_from_excluded_suite(d):
    for f in audit.pglob(os.path.join(d, 'dev/tests/test_*.m')):
        s = open(f, encoding='utf-8').read()
        if 'lsqvce.engine.resolveNc(' in s:
            open(f, 'w', encoding='utf-8').write(
                s.replace('lsqvce.engine.resolveNc(', 'zzzGone('))
    append(os.path.join(d, 'dev/tests/test_speed.m'),
           '\nfunction testZzzSpeedOnly(tc)\n% categories: contract\n'
           '    lsqvce.engine.resolveNc({struct(\'nc\', 1)});\nend\n')


case('16 a unit credited only by a comment fires', audit.check_test_categories, 'cat',
     credit_by_comment_only)
case('16 a unit credited only by the excluded benchmark fires', audit.check_test_categories,
     'cat', credit_from_excluded_suite)
case('16 a stale exemption row fires', audit.check_test_categories, 'cat',
     lambda d: append(os.path.join(d, 'dev/tests/category_exemptions.csv'),
                      'kernel.robustWeights,contract,zzz selftest\n'))
case('16 a contract exemption for a validating unit fires', audit.check_test_categories,
     'cat', lambda d: append(os.path.join(d, 'dev/tests/category_exemptions.csv'),
                             'engine.runMeta,contract,zzz selftest\n'))
case('16 a duplicated exemption row fires', audit.check_test_categories, 'cat',
     lambda d: append(os.path.join(d, 'dev/tests/category_exemptions.csv'),
                      'engine.QRSystem,speed,zzz twice\n'))

# ---- R7-32: a fault case for every check that had none --------------------
PKGF = '+lsqvce/+util/getfielddef.m'
ZZZ_SUITE = ('function tests = test_zzz_selftest\n'
             '    tests = functiontests(localfunctions);\nend\n'
             'function testZzz(tc)\n% categories: contract\n'
             '    tc.verifyTrue(true);\nend\n')
case('1 an unregistered suite fires', audit.check_test_registration, 'tests',
     lambda d: append(os.path.join(d, 'dev/tests/test_zzz_selftest.m'), ZZZ_SUITE))
case('4 a unit without its html page fires', audit.check_doc_coverage, 'docs',
     lambda d: os.remove(os.path.join(d, 'html/lsqvce.fit.html')))
case('5 a dead html link fires', audit.check_links, 'links',
     lambda d: append(os.path.join(d, 'html/index.html'), '<a href="zzz_nowhere.html">x</a>\n'))
case('6 a .js file fires', audit.check_no_js, 'js',
     lambda d: append(os.path.join(d, 'html/zzz.js'), 'var x = 1;\n'))
case('8 a hardcoded banner fires', audit.check_banners, 'banner',
     lambda d: append(os.path.join(d, PKGF),
                      "\nfunction zzzBanner()\n    lsqvce.util.log.msg('[zzz] hi\\n');\nend\n"))
case('10b a second makeKeyFcn fires', audit.check_dedup_lock, 'dedup-lock',
     lambda d: append(os.path.join(d, TFILE), '\nfunction k = makeKeyFcn(x)\n    k = x;\nend\n'))
case('11 retired stage vocabulary in output fires', audit.check_deprecated_words, 'words',
     lambda d: append(os.path.join(d, PKGF),
                      "\nfunction zzzWords()\n    lsqvce.util.log.msg('OUTER PASS 1\\n');\nend\n"))
case('12 an addpath to a missing folder fires', audit.check_addpath_targets, 'addpath',
     lambda d: append(os.path.join(d, TFILE),
                      "\nfunction zzzPath(root)\n    addpath(fullfile(root, 'zzz_nope'));\nend\n"))
case('13 a removed name in an html example fires', audit.check_html_examples, 'html-example',
     lambda d: append(os.path.join(d, 'html/index.html'), '<pre>runtests(\'x\')</pre>\n'))
case('14 a cited id nothing raises fires', audit.check_error_ids, 'error-id',
     lambda d: append(os.path.join(d, PKGF), '\n% lsqvce:zz:ghost is cited here\n'))
case('14 an id raised only by a TEST does not count', audit.check_error_ids, 'error-id',
     lambda d: (append(os.path.join(d, PKGF), '\n% lsqvce:zz:ghost is cited here\n'),
                append(os.path.join(d, TFILE),
                       "\nfunction zzzGhost()\n    error('lsqvce:zz:ghost', 'x');\nend\n")))
case('17 a wrong test-point claim fires', audit.check_readme_sync, 'readme',
     lambda d: append(os.path.join(d, TFILE),
                      '\nfunction testZzzExtra(tc)\n% categories: contract\n'
                      '    tc.verifyTrue(true);\nend\n'))
case('17 a block-commented test is not counted', audit.check_readme_sync, 'readme',
     lambda d: append(os.path.join(d, TFILE),
                      '\n%{\nfunction testZzzHidden(tc)\n% categories: contract\n'
                      '    tc.verifyTrue(true);\nend\n%}\n'),
     expect_fire=False)
case('18 a tier row naming nothing fires', audit.check_tiers_and_case, 'tiers',
     lambda d: append(os.path.join(d, 'dev/docs_src/tiers.csv'), 'zzzGhost,USER,selftest\n'))
case('22 a page that drifted from its help fires', audit.check_pages_match_help, 'docdrift',
     lambda d: append(os.path.join(d, 'html/lsqvce.fit.html'), '<p>zzz drift</p>\n'))
def manual_knobs_only():
    # check 23's rebuild half needs the CI runner's pandoc; this case tests
    # its OTHER half (every knob and group field described), so a box without
    # pandoc (the Windows baseline in CLAUDE.md) still runs the self-test clean
    audit.check_manual_built_and_complete()
    audit.FINDINGS[:] = [f for f in audit.FINDINGS if 'could not be rebuilt' not in f[1]]


case('23 a group field missing from the API reference fires', manual_knobs_only,
     'manual', lambda d: edit(os.path.join(d, 'dev/docs_src/api/04_options.md'),
                              '`SkipThreshold`', 'SkipThreshold'))


# ---- 27 groups (R7-22) and 30 (R7-28): the new instruments ---------------
def add_dead_group_field(d):
    edit(os.path.join(d, '+lsqvce/Options.m'),
         "'SkipThreshold', 0.2);", "'SkipThreshold', 0.2, 'ZzzDead', 1);")


case('27 a group field nothing reads fires', audit.check_options_are_consumed, 'options',
     add_dead_group_field)
def options_with_injected_debt():
    # the fault is a debt entry for a field the package READS (WP-116 emptied
    # the real list); the marker file in the scratch tree switches it on
    saved = set(audit.DEAD_GROUP_DEBT)
    if os.path.exists('zzz_debt_marker'):
        audit.DEAD_GROUP_DEBT.add(('Noise', 'Grid'))
    try:
        audit.check_options_are_consumed()
    finally:
        audit.DEAD_GROUP_DEBT.clear()
        audit.DEAD_GROUP_DEBT.update(saved)


case('27 a debt entry for a field that is read fires', options_with_injected_debt, 'options',
     lambda d: append(os.path.join(d, 'zzz_debt_marker'), 'x'))


def add_uncatalogued_untested_id(d):
    append(os.path.join(d, PKGF),
           "\nfunction zzzNew()\n    error('lsqvce:getfielddef:zzzNew', 'x');\nend\n")
    append(os.path.join(d, 'dev/docs_src/error_ids.csv'),
           'lsqvce:getfielddef:zzzNew,util.getfielddef,error,v6.1\n')


case('30 a new catalogued id no test fires fires', audit.check_ids_fired_by_tests, 'id-fired',
     add_uncatalogued_untested_id)
case('30 a baseline row for an id a test fires fires', audit.check_ids_fired_by_tests,
     'id-fired', lambda d: append(os.path.join(d, 'dev/docs_src/unfired_ids_baseline.csv'),
                                  'lsqvce:fit:badPrior,fit\n'))

# ---- 28 (R7-25): wrapped help lines and runtime message bodies -----------
ZZZ_UNIT = '+lsqvce/+util/zzzSelftest.m'
case('28 a claim wrapped across two help lines fires', audit.check_superseded_claims,
     'superseded', lambda d: append(os.path.join(d, ZZZ_UNIT),
                                    'function zzzSelftest()\n%ZZZSELFTEST  x.\n'
                                    '%   fewer live rows than parameters + 1 (beyond the\n'
                                    '%   breakdown point; selftest)\nend\n'))
case('28 a claim inside a runtime message fires', audit.check_superseded_claims,
     'superseded', lambda d: append(os.path.join(d, ZZZ_UNIT),
                                    "function zzzSelftest()\n%ZZZSELFTEST  x.\n"
                                    "    error('lsqvce:zzz:x', ['left too few rows: beyond the ' ...\n"
                                    "          'breakdown point of any M-estimator']);\nend\n"))
case('28 the same words split by a real line break in prose are silent',
     audit.check_superseded_claims, 'superseded',
     lambda d: append(os.path.join(d, ZZZ_UNIT),
                      'function zzzSelftest()\n%ZZZSELFTEST  x.\n'
                      '%   beyond the scope here.\n%   The breakdown is elsewhere.\nend\n'),
     expect_fire=False)

# ---- 31 (R7-35): an accepted Name-Value option its help never names --------
ZZZ_NV = '+lsqvce/+engine/zzzSelftestNV.m'
case('31 an undocumented accepted option fires', audit.check_nv_options_documented, 'nvdoc',
     lambda d: append(os.path.join(d, ZZZ_NV),
                      "function zzzSelftestNV(opt)\n%ZZZSELFTESTNV  x.\n"
                      "%   'Alpha'  the documented one\n%\n%   See also LSQVCE.FIT\n"
                      "arguments\n    opt.Alpha (1,1) double = 1\n"
                      "    opt.Hidden (1,1) double = 2\nend\nend\n"))
case('31 a documented option is silent', audit.check_nv_options_documented, 'nvdoc',
     lambda d: append(os.path.join(d, ZZZ_NV),
                      "function zzzSelftestNV(opt)\n%ZZZSELFTESTNV  x.\n"
                      "%   'Alpha'  the documented one\n%\n%   See also LSQVCE.FIT\n"
                      "arguments\n    opt.Alpha (1,1) double = 1\nend\nend\n"),
     expect_fire=False)

# ---- 32 (R7-34): a reference to a unit that does not exist ----------------
case('32 a See-also entry naming a deleted unit fires', audit.check_unit_references_resolve,
     'unitref', lambda d: append(os.path.join(d, ZZZ_UNIT),
                                 'function zzzSelftest()\n%ZZZSELFTEST  x.\n%\n'
                                 '%   See also LSQVCE.ENGINE.ADJUSTVCE\nend\n'))
case('32 a runtime remedy naming a deleted unit fires', audit.check_unit_references_resolve,
     'unitref', lambda d: append(os.path.join(d, ZZZ_UNIT),
                                 "function zzzSelftest()\n%ZZZSELFTEST  x.\n"
                                 "    error('lsqvce:zzz:x', 'use lsqvce.engine.robust instead');\nend\n"))
case('32 a deleted unit inside a dated stamp is history, silent',
     audit.check_unit_references_resolve, 'unitref',
     lambda d: append(os.path.join(d, ZZZ_UNIT),
                      'function zzzSelftest()\n%ZZZSELFTEST  x.\n%\n%   See also LSQVCE.FIT\n%\n'
                      '%   Updated 2026-10-10, lsqvce.engine.adjust retired.\nend\n'),
     expect_fire=False)


# ---- 33 (R7-44): long functions / long lines are recorded ceilings ---------
def add_long_function(d):
    append(os.path.join(d, '+lsqvce/+util/zzzLong.m'),
           'function zzzLong()\n%ZZZLONG  x.\n' + '    x = 1;\n' * 420 + 'end\n')
def grow_recorded_function(d):
    # WP-86 part 4: the recorded functions live in verdict.m; three code lines in
    # front of snoop's header land in diagnose's measured region (465 recorded)
    p = os.path.join(d, '+lsqvce/+engine/verdict.m')
    s = open(p, encoding='utf-8').read()
    s = s.replace('\nfunction [sol, D, REJ] = snoop(', '\n    x = 1;\n' * 3 + 'function [sol, D, REJ] = snoop(', 1)
    open(p, 'w', encoding='utf-8').write(s)
case('33 an unrecorded long function fires', audit.check_length_ratchets, 'length', add_long_function)
case('33 a recorded long function that grows fires', audit.check_length_ratchets, 'length',
     grow_recorded_function)
case('33 a stale long-function row fires', audit.check_length_ratchets, 'length',
     lambda d: append(os.path.join(d, 'dev/docs_src/long_functions.csv'),
                      '+lsqvce/fit.m,fit,999,"no longer long"\n'))
case('33 a new line over 120 characters fires', audit.check_length_ratchets, 'length',
     lambda d: append(os.path.join(d, ZZZ_UNIT),
                      'function zzzSelftest()\n%ZZZSELFTEST  x.\n    y = ' + '1 + ' * 40 + '1;\nend\n'))


# ---- 34 (WP-138): the ACCURACY ratchet ----------------------------------------
case('34 a new unit without an ACCURACY block fires', audit.check_accuracy_ratchet, 'accuracy',
     lambda d: append(os.path.join(d, ZZZ_UNIT),
                      'function zzzSelftest()\n%ZZZSELFTEST  x.\n%\n%   See also LSQVCE.FIT\nend\n'))
case('34 a new unit that states its ACCURACY is silent', audit.check_accuracy_ratchet, 'accuracy',
     lambda d: append(os.path.join(d, ZZZ_UNIT),
                      'function zzzSelftest()\n%ZZZSELFTEST  x.\n%\n%   ACCURACY\n'
                      '%     Not applicable: no arithmetic.\nend\n'),
     expect_fire=False)
case('34 a frozen row whose unit states its ACCURACY is stale and fires',
     audit.check_accuracy_ratchet, 'accuracy',
     lambda d: append(os.path.join(d, 'dev/docs_src/accuracy_baseline.csv'), 'engine.verdict,v6.24.0\n'))


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
