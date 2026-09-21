#!/usr/bin/env python3
"""Audit the lsqvce tree for internal consistency.

WHY THIS EXISTS
    The HANDOVER can drift out of sync with the code (it did, on 2026-07-18:
    six tests were calling a removed signature, twelve functions had no doc
    page, and several help blocks still described deleted aliases). Trusting a
    prose handover is not a verification method. This script re-derives the
    facts from the tree, so the handover can be CHECKED instead of believed.

    It is READ-ONLY: it never edits anything, it only reports.

WHAT IT CHECKS
    1  every test file on disk is registered in tests/run_all_lsqvce.m
       (the suite list is hardcoded, so a new file is invisible until added)
    2  removed names do not appear in executable code
    3  public entry-point call sites match the current signatures
    4  every +lsqvce function has an html/ page and a Contents.m entry
    5  every html link and helptoc target resolves
    6  the doc set stays JavaScript-free (hard environment requirement)
    7  function/end balance, with MATLAB-aware tokenising
    ...
    18 surface census (dev/docs_src/tiers.csv complete both ways) and no
       case-only file-name collisions (Windows working copies)
    19 catch-block ratchet: every empty or bare `catch` (no identifier
       inspected, everything swallowed) is listed in
       dev/docs_src/catch_baseline.csv; a NEW one is a finding, a stale row
       is a finding (R5 M0; the doctrine target is an empty baseline, M1)
    20 error-identifier catalogue: every lsqvce:* id the package raises is a
       row of dev/docs_src/error_ids.csv and vice versa (closed both ways);
       ids raised but not documented in their own unit's help are REPORTED
       (the M1 debt number)
    21 layering ratchet: tiers.csv gives every public name a tier
       (USER > ENGINE > KERNEL > UTIL); calls go DOWNWARD only. Every upward
       edge that exists today is listed in dev/docs_src/layering_baseline.csv;
       a NEW upward edge is a finding, a stale row is a finding
    22 html pages match the help blocks they were generated from
    23 the Manual is a FRESH build of dev/docs_src/manual/*.md (byte equality;
       needs the same pandoc as the CI runner)
    24 help EXAMPLE blocks call units that exist
    25 every unit page is reachable from html/index.html and helptoc.xml
    26 markdown pointers resolve: repository-relative links and backticked
       paths in the documents that TELL A READER WHERE TO GO (R6 WP-99;
       dev/HANDOVER.md and dev/HISTORY_v1.md are records and are excluded --
       naming a deleted file is half their job)
    29 the PDF Manual's title carries the Contents.m version (check 23
       rebuilds the HTML Manual and compares bytes; the PDF cannot be
       compared that way, so the ONE fact that must not drift is read)
    27 every public Options property has a consumer: engineArgs, or somewhere
       in +lsqvce outside Options.m (R6 WP-101; MinPasses and LogFile shipped
       dead for a whole review)
    THE AUTHORITATIVE LIST is the order of calls in main(); this summary is a
       reading aid and was already four checks behind when 26 was added.
    Self-test: dev/docs_src/audit_selftest.py injects one fault per new
       check into a scratch copy and asserts the finding fires (CI runs it).

USAGE
    python3 dev/docs_src/audit_tree.py            # from the toolbox root
    echo $?                                   # 0 = clean, 1 = findings

Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""

import binascii
import glob
import os
import re
import sys
import zlib
import xml.etree.ElementTree as ET


def pglob(pattern, **kw):
    """glob.glob with POSIX separators on every platform (R6 WP-93, E-12).
    Every check derives unit names, tiers and index keys by string work on
    '/'-separated paths; on Windows glob returns backslashes, which turned
    check 16 into ~590 false gaps and check 21 into stale-row findings. One
    helper, one convention: the audit reads the same tree the same way
    everywhere. Forward slashes open files on Windows too."""
    return [p.replace('\\', '/') for p in glob.glob(pattern, **kw)]

FINDINGS = []


def note(check, msg):
    FINDINGS.append((check, msg))
    print('  FINDING [%s] %s' % (check, msg))


# --------------------------------------------------------------------------
def strip_matlab(line):
    """Blank out comments and string literals, keeping transposes intact."""
    out = []
    i, n, instr = 0, len(line), None
    while i < n:
        c = line[i]
        if instr:
            out.append(' ')
            if c == instr:
                instr = None
            i += 1
            continue
        if c == '%':
            break
        if line.startswith('...', i):
            break
        if c == '"':
            instr = '"'
            out.append(' ')
            i += 1
            continue
        if c == "'":
            j = len(out) - 1
            while j >= 0 and out[j] == ' ':
                j -= 1
            prev = out[j] if j >= 0 else ''
            if re.match(r"[\w\)\]\}\.']", prev):
                out.append("'")            # transpose
            else:
                instr = "'"
                out.append(' ')
            i += 1
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def code_lines(path):
    """Yield (lineno, code-only text) for a .m file."""
    with open(path, encoding='utf-8', errors='replace') as fh:
        for k, raw in enumerate(fh, 1):
            s = strip_matlab(raw)
            if s.strip():
                yield k, s


# --------------------------------------------------------------------------
def check_test_registration():
    print('\n1. test-file registration')
    runner = os.path.join('dev', 'tests', 'run_all_lsqvce.m')
    if not os.path.exists(runner):
        note('tests', 'dev/tests/run_all_lsqvce.m is missing')
        return
    reg = set(re.findall(r"'(test_\w+)\.m'", open(runner, encoding='utf-8').read()))
    disk = {os.path.basename(f)[:-2]
            for f in pglob('dev/tests/*.m') + pglob('dev/tests/accum/*.m')
            if os.path.basename(f).startswith('test_')}
    disk.discard('test_speed')                 # a benchmark, run manually
    for t in sorted(disk - reg):
        note('tests', '%s is on disk but NOT in run_all_lsqvce' % t)
    for t in sorted(reg - disk):
        note('tests', '%s is registered but NOT on disk' % t)
    print('   %d suites registered, %d on disk' % (len(reg), len(disk)))


def check_removed_names():
    print('\n2. removed names in executable code')
    removed = ['SeparableProblem', 'covModelFactory', 'fromAccum', 'psdIterate',
               'descaledResidualAcv', 'maxLagCorr', 'selectAROrder',
               'vceGroupsBy', 'coalesceOpt']
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            for k, s in code_lines(f):
                for nm in removed:
                    if re.search(r'\b%s\b' % nm, s):
                        note('names', '%s:%d uses removed name %r' % (f, k, nm))
                if re.search(r'\blsqvce\.recover\s*\(', s):
                    note('names', '%s:%d calls removed lsqvce.recover' % (f, k))

    # ALSO CHECK HELP TEXT. code_lines() strips comments, so the loop above
    # cannot see a `See also LSQVCE.SEPARABLEPROBLEM` -- and MATLAB renders those
    # as clickable links that then fail. Fifteen such references survived the
    # rename undetected, including in the two demo files a user is told to copy.
    # This is the third defect traced to "the checker could not see comments";
    # read the RAW source here.
    hist = ('REMOVED', 'Renamed 2026', 'deprecated', 'Was SeparableProblem',
            'was SeparableProblem')
    dead = re.compile(r'SeparableProblem|SEPARABLEPROBLEM|LSQVCE\.RECOVER'
                      r'|lsqvce\.recover\b', re.I)
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            for k, raw in enumerate(open(f, encoding='utf-8',
                                         errors='replace').read().split('\n'), 1):
                if not raw.lstrip().startswith('%'):
                    continue                    # code is covered above
                if not dead.search(raw):
                    continue
                if any(h in raw for h in hist):
                    continue                    # a deliberate history note
                note('names', '%s:%d help text still names a removed function'
                     % (f, k))


def check_signatures():
    print('\n3. entry-point signatures')
    # nc was removed as a leading argument on 2026-07-18.
    pat = re.compile(r'\badjust(?:Block|VCE)\s*\(\s*(?:nc\b|\d+\s*,)')
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            for k, s in code_lines(f):
                if pat.search(s):
                    note('signature',
                         '%s:%d passes nc as a leading argument (removed)'
                         % (f, k))

    # Deriving nc happens EARLY unless the entry point validates the SHAPE of
    # `systems` first, so a validation test that feeds a deliberately malformed
    # `systems` can be pre-empted by resolveNc and get the generic
    # lsqvce:resolveNc:missing instead of the specific id it asserts. This is
    # how testRecoverErrors broke: verifyError(@() adjustBlock({1},1:3),
    # 'lsqvce:adjustBlock:badCell') hit resolveNc first.
    #
    # Only warn when adjustBlock ACTUALLY still derives nc before validating,
    # so the check reports the defect and goes quiet once it is fixed. A check
    # that cries wolf on a healthy tree gets ignored (HANDOVER §4d).
    src = ''
    if os.path.exists('+lsqvce/adjustBlock.m'):
        src = '\n'.join(s for _, s in code_lines('+lsqvce/adjustBlock.m'))
    i_nc = src.find('resolveNc')
    i_norm = src.find('normalizeSystems(')
    derives_first = (0 <= i_nc < i_norm) if (i_nc >= 0 and i_norm >= 0) else False

    if derives_first:
        # NOTE: reads the RAW source, not code_lines(). code_lines() blanks
        # string literals (correctly, for balance counting), which would erase
        # the very error id being matched. The leading segment is [^;]*? so a
        # match cannot span past a statement end into an unrelated assertion.
        idpat = re.compile(
            r"verifyError\s*\([^;]*?adjust(?:Block|VCE)\s*\(([^;]*?)\)\s*,\s*"
            r"\.{0,3}\s*'(lsqvce:[\w:]+)'", re.S)
        for f in pglob('dev/tests/**/*.m', recursive=True):
            raw = open(f, encoding='utf-8', errors='replace').read()
            for m in idpat.finditer(raw):
                arg, want = m.group(1), m.group(2)
                if want.endswith(':empty') or 'resolveNc' in want:
                    continue                   # resolveNc's own ids are intended
                if not re.search(r"\bnc\b", arg):
                    note('signature',
                         '%s: verifyError expects %s, but adjustBlock derives '
                         'nc BEFORE normalizeSystems validates the shape, so '
                         'resolveNc:missing fires first' % (f, want))


def check_doc_coverage():
    print('\n4. doc coverage')
    # R4-D1 (WP-64a): the package has sub-packages +engine / +kernel / +util;
    # a name is listed in ITS package's Contents.m and its page is
    # html/lsqvce.<pkg>.<name>.html.
    funcs = []
    for pkg, sub in (('', ''), ('engine.', '+engine/'), ('kernel.', '+kernel/'),
                     ('util.', '+util/')):
        contents = ''
        cpath = '+lsqvce/%sContents.m' % sub
        if os.path.exists(cpath):
            contents = open(cpath, encoding='utf-8').read()
        elif sub:
            note('docs', '%s missing' % cpath)
        for f in sorted(pglob('+lsqvce/%s*.m' % sub)):
            b = os.path.basename(f)[:-2]
            if b != 'Contents':
                funcs.append((pkg + b, 'html/lsqvce.%s%s.html' % (pkg, b), contents, cpath))
    for name, page, contents, cpath in funcs:
        if not os.path.exists(page):
            note('docs', 'no html page for %s (doc lsqvce.%s fails)' % (name, name))
        short = name.split('.')[-1]
        if not re.search(r'\b%s\b' % re.escape(short), contents):
            note('docs', '%s is not listed in %s' % (name, cpath))
    print('   %d functions checked' % len(funcs))

    for page in sorted(pglob('html/lsqvce.*.html')):
        n = os.path.basename(page)[len('lsqvce.'):-len('.html')]
        m = '+lsqvce/%s.m' % n
        for pkg, sub in (('util.', '+util/'), ('engine.', '+engine/'), ('kernel.', '+kernel/')):
            if n.startswith(pkg):
                m = '+lsqvce/%s%s.m' % (sub, n[len(pkg):])
        if not os.path.exists(m):
            # WP-41 decision (user-approved by PR merge): redirect stubs for
            # REMOVED v1 names are KEPT for old bookmarks -- but only under a
            # formal rule, enforced here: the page must (a) redirect a name
            # that appears in the removed/renamed registry below, (b) carry
            # "(removed)" in its <title>, and (c) meta-refresh to a live
            # page. Anything else without a backing .m is a FINDING.
            stubOK = {'lsqvce.recover.html', 'lsqvce.SeparableProblem.html'}
            pageBase = os.path.basename(page)
            pth = page if os.path.exists(page) else os.path.join('html', page)
            raw = open(pth, encoding='utf-8', errors='replace').read()
            isStub = ('(removed)' in raw.split('</title>')[0]
                      and 'http-equiv="refresh"' in raw)
            if pageBase in stubOK and isStub:
                print('   OK stub: %s (removed-name redirect, kept by the '
                      'WP-41 decision)' % page)
            else:
                note('doc-coverage', '%s has no backing .m and is not an '
                     'approved removed-name redirect stub' % page)


def check_links():
    print('\n5. link integrity')
    bad = 0
    for f in pglob('html/*.html'):
        s = open(f, encoding='utf-8', errors='replace').read()
        for m in re.finditer(r'href="([^"#?]+)"', s):
            t = m.group(1)
            if t.startswith('http') or t.startswith('mailto'):
                continue
            if not os.path.exists(os.path.join('html', t)):
                note('links', '%s -> %s does not exist' % (f, t))
                bad += 1
    if os.path.exists('helptoc.xml'):
        try:
            for it in ET.parse('helptoc.xml').getroot().iter('tocitem'):
                t = it.get('target')
                if t and not os.path.exists(t):
                    note('links', 'helptoc target missing: %s' % t)
                    bad += 1
        except ET.ParseError as e:
            note('links', 'helptoc.xml is not well-formed: %s' % e)
    if not bad:
        print('   all links resolve')


def check_no_js():
    print('\n6. JavaScript-free doc set')
    n = 0
    for f in pglob('html/**/*.html', recursive=True):
        s = open(f, encoding='utf-8', errors='replace').read()
        if '<script' in s.lower():
            note('js', '%s contains a <script> tag' % f)
            n += 1
    for f in pglob('**/*.js', recursive=True):
        note('js', 'JavaScript file present: %s' % f)
        n += 1
    if not n:
        print('   no <script> tags, no .js files')


def check_balance():
    print('\n7. function/end balance')
    OPEN = {'function', 'if', 'for', 'while', 'switch', 'try', 'parfor',
            'classdef', 'methods', 'properties', 'arguments', 'enumeration',
            'events', 'spmd'}
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            depth = bdepth = 0
            for _, s in code_lines(f):
                # A `case` label is DATA, not code: `case 'parfor'` and
                # `case 'spmd'` contain block-opener keywords as string
                # literals. The quote scanner reads the leading quote after
                # `case ` as a transpose (the preceding token is a word), so
                # the literal survives stripping and would be miscounted.
                # `case` opens no block, so dropping the label is safe.
                s = re.sub(r'^\s*case\b.*$', '', s)
                for m in re.finditer(r"[A-Za-z_]\w*|\(|\)|\[|\]|\{|\}", s):
                    t = m.group(0)
                    if t in '([{':
                        bdepth += 1
                    elif t in ')]}':
                        bdepth -= 1
                    elif bdepth > 0:
                        continue            # inside brackets 'end' is an index
                    elif t in OPEN:
                        depth += 1
                    elif t == 'end':
                        depth -= 1
            if depth != 0:
                note('balance', '%s is unbalanced by %+d' % (f, depth))
    print('   done')


def check_banners():
    print('\n8. banner vocabulary is centralised')
    # verbosity.banner exists so every tag renders at ONE width. That only holds
    # if nothing writes its own literal: the 2026-07-18 MATLAB log showed
    # '[resid ]' 499 times and '[resid   ]' 63 times -- the same tag at two
    # widths in a single run, because only adjustVCE and BlockAccumulator's loop
    # had been converted. This flags any log.msg format string that still embeds
    # a bracketed word tag.
    #
    # NOT flagged: units ('[Hz]', '[s]', '[%.2f s]'), report.m's severity
    # markers ('[ERROR]', '[WARN ]') which are a deliberately narrow column in a
    # text table, and verbosity.m itself, whose help shows examples.
    allow_files = {'+lsqvce/+util/verbosity.m', '+lsqvce/report.m'}
    unit_tags = {'Hz', 's', 't', 'm', 'k'}
    pat = re.compile(r"log\.(?:msg|warn)\(\s*\[?\s*'[^']*?\[([A-Za-z]+) *\]")
    n = 0
    for f in pglob('+lsqvce/**/*.m', recursive=True):
        if f.replace('\\', '/') in allow_files:
            continue
        # RAW source, not code_lines(): code_lines blanks string literals
        # (correctly, for balance counting), which would erase the very banner
        # text being matched. This exact trap already cost one silent miss.
        for k, raw in enumerate(open(f, encoding='utf-8',
                                     errors='replace').read().split('\n'), 1):
            s = raw.split('%')[0] if raw.lstrip().startswith('%') else raw
            for m in pat.finditer(s):
                if m.group(1) in unit_tags:
                    continue
                note('banner', '%s:%d hardcodes the banner [%s] instead of '
                     'lsqvce.util.verbosity.banner' % (f, k, m.group(1)))
                n += 1
    if not n:
        print('   no hardcoded banner literals')


def check_spec_fields():
    print('\n9. block-spec field list agrees everywhere')
    # The spec struct is written out in FOUR places: a repmat template and a
    # struct(...) literal in each of adjustBlock and adjustVCE. MATLAB rejects
    # `specs(i) = struct(...)` when the field sets differ, and it does so at RUN
    # time -- nothing catches it until the loop executes. Adding 'name' on
    # 2026-07-19 updated two of the four and broke demo_messaging.
    pat = re.compile(r"struct\(\s*((?:'[A-Za-z_]\w*'\s*,[^;]*?))\)", re.S)
    sets = {}
    for f in ('+lsqvce/adjustBlock.m', '+lsqvce/adjustVCE.m'):
        if not os.path.exists(f):
            continue
        raw = open(f, encoding='utf-8', errors='replace').read()
        for m in re.finditer(r"^[^%\n]*\b(?:tmpl|specs\(i\)|s)\s*=\s*struct\(",
                             raw, re.M):
            seg = raw[m.end() - 7:]
            mm = pat.match(seg)
            if not mm:
                continue
            names = tuple(re.findall(r"'([A-Za-z_]\w*)'\s*,", mm.group(1)))
            # keep only the ones that look like a block spec
            if 'sysID' in names and 'dayList' in names:
                sets.setdefault(tuple(sorted(names)), []).append(
                    '%s:%d' % (f, raw[:m.start()].count('\n') + 1))
    if len(sets) > 1:
        for names, where in sets.items():
            note('spec', 'block-spec fields %s at %s'
                 % (', '.join(names), '; '.join(where)))
    else:
        n = len(next(iter(sets.values()))) if sets else 0
        print('   %d definition(s), one field list' % n)


def check_dedup_lock():
    print('\n10b. dedup lock: ONE makeKeyFcn / ONE arcBatchWhiten (WP-38, F-08)')
    # These two skeletons existed as verbatim copies in redundancy.m and
    # streamContribution.m until WP-38; a silent re-copy would let the
    # grouping keys or the batch/whitening semantics drift apart again.
    # The lock: each name is DEFINED exactly once, and only under
    # +lsqvce/private/.
    pat = re.compile(r'^\s*function\s+(?:(?:\[[^\]]*\]|[\w~]+)\s*=\s*)?'
                     r'(makeKeyFcn|arcBatchWhiten)\s*\(', re.M)
    hits = {}
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            raw = open(f, encoding='utf-8', errors='replace').read()
            code = '\n'.join(l for l in raw.split('\n')
                             if not l.lstrip().startswith('%'))
            for m in pat.finditer(code):
                hits.setdefault(m.group(1), []).append(f)
    ok = True
    for name in ('makeKeyFcn', 'arcBatchWhiten'):
        where = hits.get(name, [])
        want = os.path.join('+lsqvce', '+engine', 'private', name + '.m')   # R4-D1: engine-private since WP-64a
        if len(where) != 1 or os.path.normpath(where[0]) != os.path.normpath(want):
            note('dedup-lock', '%s must be defined EXACTLY once, in %s (found: %s)'
                 % (name, want, ', '.join(where) or 'nowhere'))
            ok = False
    if ok:
        print('    OK - both defined once, under +lsqvce/+engine/private/')

def check_duplicate_functions():
    print('\n10. no duplicated function names within a file')
    # MATLAB refuses to load a file that defines the same function twice, so the
    # whole suite file silently contributes NOTHING -- and run_all_lsqvce
    # reported "0 failed" while 41 tests never ran. A duplicate is easy to
    # create when editing a long test file by inserting a block, which is
    # exactly how it happened on 2026-07-19.
    # The output-argument list and its '=' are ONE optional unit. Written as
    # two separate optionals, the name group happily matched the LAST letter of
    # the function name ('testFoo(' -> 'o'), so every file looked duplicated.
    pat = re.compile(r'^\s*function\s+(?:(?:\[[^\]]*\]|[\w~]+)\s*=\s*)?'
                     r'([A-Za-z]\w*)\s*\(')
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            raw = open(f, encoding='utf-8', errors='replace').read()
            names = []
            for line in raw.split('\n'):
                if line.lstrip().startswith('%'):
                    continue
                m = pat.match(line)
                if m:
                    names.append(m.group(1))
            seen, dup = set(), []
            for n in names:
                if n in seen and n not in dup:
                    dup.append(n)
                seen.add(n)
            for n in dup:
                note('duplicate', '%s defines %s() more than once - MATLAB '
                     'cannot load this file at all' % (f, n))


def check_deprecated_words():
    print('\n11. no deprecated stage vocabulary in user-visible text')
    # SCALE / SHAPE / SOLVE / GATE were INVENTED stage names, rejected by the
    # user as non-geodetic and replaced by plain words on 2026-07-18. They kept
    # reappearing in output because a rename that fixes the banner does not fix
    # the sentence next to it: the 'OUTER PASS ... SCALE ... SOLVE ... SHAPE ...
    # GATE' sub-line survived two rounds, and the notConverged WARNING still
    # said 'both gates' a round after that.
    #
    # Only strings that a USER can see are checked: log.msg, log.warn, error,
    # warning, title, xlabel, ylabel, fprintf. Ordinary comments are prose and
    # may legitimately say 'the shape of the noise'.
    emit = re.compile(r'\b(?:log\.msg|log\.warn|error|warning|title|xlabel|'
                      r'ylabel|fprintf)\s*\(')
    dead = re.compile(r"\b(SCALE|SHAPE|GATE|SOLVE|OUTER PASS|both gates)\b")
    allow = ('verbosity.m',)          # its help QUOTES the retired form
    for d in ('+lsqvce', 'tests', 'examples'):
        for f in pglob(d + '/**/*.m', recursive=True):
            if os.path.basename(f) in allow:
                continue
            raw = open(f, encoding='utf-8', errors='replace').read().split('\n')
            # JOIN CONTINUATIONS FIRST. The deprecated words lived on the SECOND
            # line of a wrapped log.msg call, so a line-by-line check that
            # required the call and the word on one line saw nothing -- which is
            # exactly how they survived three rounds of searching.
            k = 0
            while k < len(raw):
                start = k
                line = raw[k]
                while line.rstrip().endswith('...') and k + 1 < len(raw):
                    k += 1
                    line = line.rstrip()[:-3] + ' ' + raw[k].lstrip()
                k += 1
                if raw[start].lstrip().startswith('%'):
                    continue          # a comment is prose, not output
                if not emit.search(line):
                    continue
                for m in dead.finditer(line):
                    before = line[:m.start()]
                    if before.count("'") % 2 == 1:   # inside a quoted literal
                        note('words', '%s:%d user-visible text still says %r'
                             % (f, start + 1, m.group(1)))


def check_addpath_targets():
    print('\n12. every addpath target exists')
    # lsqvce_setup.m added THREE folders that had been deleted -- tests/accum
    # (removed by the test flattening, HANDOVER 1.7b) and examples/accum +
    # examples/gravity (removed by the examples consolidation, 1.3b), both on
    # 2026-07-18. Every user running setup got three
    # MATLAB:addpath:nameNonexistentOrNotADirectory warnings, and it went
    # unnoticed until 2026-07-19 because nothing re-derived the path list from
    # the tree. Same class as HANDOVER 1.0 defect 4: a removal recorded as done
    # with one surface left behind.
    #
    # Only LITERAL fullfile(root, ...) targets can be checked statically; a
    # computed path is skipped rather than guessed at.
    call = re.compile(r"addpath\s*\(\s*fullfile\s*\(\s*root\s*,\s*([^)]*)\)")
    for f in pglob('*.m') + pglob('dev/tests/*.m'):
        raw = open(f, encoding='utf-8', errors='replace').read().split('\n')
        for i, line in enumerate(raw):
            if line.lstrip().startswith('%'):
                continue                       # a comment is not a call
            m = call.search(line)
            if not m:
                continue
            parts = re.findall(r"'([^']*)'", m.group(1))
            if len(parts) != len(m.group(1).split(',')):
                continue                       # a computed segment -- skip
            target = os.path.join(*parts) if parts else ''
            if target and not os.path.isdir(target):
                note('addpath', '%s:%d addpath target %r does not exist'
                     % (f, i + 1, target))


def check_html_examples():
    print('\n13. shipped HTML examples use live names')
    # TWICE in two sessions, shipped clickable docs told users to call a
    # signature that errors: 12 examples kept the removed leading `nc` argument
    # (2026-07-19, HANDOVER 1.1d) and 14 kept the pre-rename fluent methods
    # .vceGroupsBy/.shape/.scale (2026-07-19, HANDOVER 1.2a). Check 3 covers
    # entry-point signatures in .m files ONLY, so both survived a clean audit.
    # A doc example is executable text: a user copies it verbatim.
    #
    # Keyed on the OLD name; the value is what it became, so the finding tells
    # the reader what to write instead of merely that something is wrong.
    removed = {
        'vceGroupsBy':  "'GroupBy' constructor group (2.0)",
        r'\.shape\(':   "'Shape' constructor group (2.0)",
        r'\.scale\(':   "'Scale' constructor group (2.0)",
        r'\.groupComponentsBy\(':  "'GroupBy', key on the constructor (R4-D2, 2.0)",
        r'\.stochasticModel\(':    "'Shape', struct(...) on the constructor (R4-D2, 2.0)",
        r'\.varianceComponents\(': "'Scale', struct(...) on the constructor (R4-D2, 2.0)",
        r'\.componentsPerEpoch\(': "'Axes', struct(...) on the constructor (R4-D2, 2.0)",
        r'adjustVCE\(nc,':    'adjustVCE(systems, days) -- nc is derived',
        r'adjustBlock\(nc,':  'adjustBlock(systems, days) -- nc is derived',
        r'adjustBlock\(\d+,': 'adjustBlock(systems, days) -- nc is derived',
        r'P\.convergence\(':  'R.convergence() -- moved to Results',
        r'runtests\(':        'run_all_lsqvce',
    }
    # Only <pre> blocks: prose may legitimately DISCUSS an old name (e.g. a
    # migration note), but an example is what gets copied.
    pre = re.compile(r'<pre\b[^>]*>(.*?)</pre>', re.S)
    for f in sorted(pglob('html/*.html')):
        body = open(f, encoding='utf-8', errors='replace').read()
        for block in pre.findall(body):
            for pat, now in removed.items():
                if re.search(pat, block):
                    note('html-example',
                         '%s example uses removed %r -- now %s'
                         % (f, pat.replace('\\', ''), now))


def check_error_ids():
    print('\n14. documented error identifiers exist in the code')
    # An identifier in help text or a doc page is a PROMISE: the user greps for
    # it, or writes a try/catch against it. If no code throws it, that promise
    # is broken silently -- nothing else in this audit would notice, because the
    # id is well-formed and the file parses.
    thrown = set()
    for f in pglob('+lsqvce/**/*.m', recursive=True) + \
             pglob('template/*.m') + pglob('dev/tests/*.m'):
        src = open(f, encoding='utf-8', errors='replace').read()
        # Every way this toolbox raises an identified condition. The log
        # wrapper (lsqvce.util.log.warn/warnOnce) is the big one: most WARNING
        # ids go through it, and a checker that only knows error() reports the
        # whole warning vocabulary as undefined -- 80 false positives, which is
        # worse than no check because it trains you to ignore the output.
        # MATLAB accepts both quote styles, and resolveNc.m uses double quotes.
        src = src.replace('("', "('").replace('",', "',")   # normalise for scanning
        for pat in (r"error\(\s*'([\w]+:[\w:]+)'",
                    r"MException\(\s*'([\w]+:[\w:]+)'",
                    r"warning\(\s*'([\w]+:[\w:]+)'",
                    r"log\.warn\w*\(\s*'([\w]+:[\w:]+)'",
                    r"warnOnce\(\s*'([\w]+:[\w:]+)'",
                    r"'([\w]+:[\w:]+)'\s*,\s*\.\.\."):
            thrown |= set(re.findall(pat, src))
    # MATLAB's own identifiers are not ours to define
    ok = lambda i: i.startswith('lsqvce')
    cited = {}
    for f in pglob('+lsqvce/**/*.m', recursive=True) + \
             pglob('html/*.html') + pglob('template/*.m'):
        for _, ln, i in [(0, n, m) for n, line in
                         enumerate(open(f, encoding='utf-8',
                                        errors='replace'), 1)
                         for m in re.findall(r'\blsqvce:[\w:]+', line)]:
            cited.setdefault(i, (f, ln))
    for i, (f, ln) in sorted(cited.items()):
        if ok(i) and i not in thrown and not i.endswith(':'):
            note('error-id', '%s:%d cites %s, which no code throws' % (f, ln, i))


def check_h1_names():
    print('\n15. H1 help lines match their filename')
    # A rename that updates identifiers but not the H1 COMMENT leaves
    # `help myFile` announcing the old name. Three files survived the
    # 2026-07-19 separable->grace rename this way (demo_grace_dataLoader,
    # demo_grace, test_adjustVCE), because check 2 guards removed FUNCTION
    # names in code, not the prose of an H1 line.
    for d in ('+lsqvce', 'tests', 'examples', 'template'):
        for f in sorted(pglob(d + '/**/*.m', recursive=True)):
            name = os.path.splitext(os.path.basename(f))[0]
            if name == 'Contents':
                continue                      # by convention names the PACKAGE
            for line in open(f, encoding='utf-8', errors='replace'):
                if line.startswith('%') and len(line) > 2:
                    w = line[1:].split()[0] if line[1:].split() else ''
                    # a package H1 is LSQVCE.NAME; compare the last component
                    if not w.replace('_', '').replace('.', '').isalnum():
                        break
                    if w.split('.')[-1].upper() != name.upper():
                        note('h1', '%s: H1 says %r' % (f, w))
                    break


def check_test_categories():
    print('\n16. test categories are declared, counted per unit, and enforced')
    # WHY (WP-43, HANDOVER par. 6.1.2 as amended 2026-08-19): the doctrine
    # demands the categories be COUNTABLE mechanically, per unit. Another
    # project's first such count found 46 of 81 functions short of at least
    # one category AFTER nine green runs. This check makes the count a
    # standing gate instead of a one-time review:
    #   (a) every test point declares '% categories: ...' on the line after
    #       its signature; only the six category names and the two
    #       cross-cutting tags are legal vocabulary;
    #   (b) the unit x category matrix is derived from the lsqvce.<name>
    #       references in each test body (mechanical, no second bookkeeping
    #       file to drift) and REPORTED every run;
    #   (c) a public unit with zero 'contract' coverage is a FINDING unless
    #       tests/category_exemptions.csv carries a reasoned row -- contract
    #       may not be exempted for a function that validates its arguments;
    #   (d) tolerance provenance: in the suites listed in TOL_SWEPT every
    #       'AbsTol'/'RelTol' literal must carry a 'tol:' comment on the
    #       same or one of the two preceding lines. Elsewhere the
    #       unannotated count is REPORTED (the debt number to drive to 0).
    # R5 M0 (WP-78, 2026-09-01): SEVEN categories -- 'metamorphic' is a
    # counted category (VALIDATION_GUIDE Part 1: invariance/equivariance
    # through the public API, the cheap one nobody writes), no longer a
    # tag; three more cross-cutting tags (layering, concurrency,
    # convergence). The matrix gains a column; gaps stay informational
    # until M7 (WP-89) drives them to 0 or to a reasoned exemption row.
    cats_ok = {'contract', 'reference', 'precision', 'speed',
               'robustness', 'vectorisation', 'metamorphic'}
    tags_ok = {'diagnostics', 'provenance', 'layering', 'concurrency',
               'convergence'}
    # WP-43 r2 (2026-08-19): the staged rollout is COMPLETE -- the debt count
    # reached 0 in the same round that annotated all 216 remaining sites, so
    # the swept-suite list retires and every suite is enforced. The debt
    # counter stays as a belt-and-braces invariant: any future unannotated
    # tolerance is a FINDING everywhere, not a reported number.
    TOL_SWEPT = sorted(pglob('dev/tests/test_*.m'))

    decl_re = re.compile(
        r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
    fn_re = re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')

    unit_cat = {}          # unit -> set of categories reaching it
    n_points = 0
    for f in sorted(pglob('dev/tests/test_*.m')):
        lines = open(f, encoding='utf-8', errors='replace').read().split('\n')
        starts = [(i, m.group(1)) for i, l in enumerate(lines)
                  if (m := fn_re.match(l))
                  and not re.match(r'function\s+tests\s*=', l)]
        for k, (i, name) in enumerate(starts):
            n_points += 1
            end = starts[k + 1][0] if k + 1 < len(starts) else len(lines)
            m = decl_re.match(lines[i + 1]) if i + 1 < len(lines) else None
            if not m:
                note('cat', '%s: %s has no "%% categories:" declaration'
                     % (f, name))
                continue
            cats = {c.strip() for c in m.group(1).split(',') if c.strip()}
            tags = {t.strip() for t in (m.group(2) or '').split(',')
                    if t.strip()}
            bad = (cats - cats_ok) | (tags - tags_ok)
            if bad or not cats:
                note('cat', '%s: %s declares unknown category/tag %s'
                     % (f, name, sorted(bad)))
                continue
            body = '\n'.join(lines[i + 1:end])
            for u in set(re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w+)', body)):
                unit_cat.setdefault(u, set()).update(cats)

    units = []
    for f in sorted(pglob('+lsqvce/**/*.m', recursive=True)):
        if '/private/' in f:
            continue
        u = (f.replace('+lsqvce/', '').replace('.m', '').replace('+util/', 'util.')
              .replace('+engine/', 'engine.').replace('+kernel/', 'kernel.'))
        if not u.endswith('Contents'):
            units.append(u)

    exempt = {}
    if os.path.isfile('dev/tests/category_exemptions.csv'):
        for row in open('dev/tests/category_exemptions.csv',
                        encoding='utf-8'):
            row = row.strip()
            if not row or row.startswith('#') or row.startswith('unit,'):
                continue
            parts = row.split(',', 2)
            if len(parts) == 3 and parts[2].strip():
                exempt[(parts[0].strip(), parts[1].strip())] = parts[2]
            else:
                note('cat', 'category_exemptions.csv row lacks a reason: %r'
                     % row)

    print('   %d test points declared; unit x category coverage '
          '(%d public units):' % (n_points, len(units)))
    hdr = '   %-22s %s' % ('unit', ' '.join('%-6s' % c[:6] for c in sorted(
        cats_ok)))
    print(hdr)
    gaps = 0
    for u in units:
        got = unit_cat.get(u, set())
        row = ' '.join('%-6s' % ('x' if c in got else '.')
                       for c in sorted(cats_ok))
        missing = sorted(cats_ok - got)
        print('   %-22s %s%s' % (u, row,
                                 ('  gap: ' + ','.join(missing))
                                 if missing else ''))
        gaps += len(missing)
        if 'contract' not in got and (u, 'contract') not in exempt:
            note('cat', 'public unit %s has NO contract coverage and no '
                 'exemption row' % u)
        # R5 M7 (WP-89): every OTHER (unit, category) gap is a finding too
        # unless an ARGUED exemption row exists -- a category is tested or
        # argued not to apply, never silently missing.
        for c in missing:
            if c != 'contract' and (u, c) not in exempt:
                note('cat', 'public unit %s has no %s coverage and no '
                     'exemption row (M7: tested or argued)' % (u, c))
    print('   category gaps in the matrix: %d (every one either covered '
          'or argued in category_exemptions.csv)' % gaps)

    # (d0) a tol: comment must never sit INSIDE a ... continuation -- a
    # standalone comment line there is a PARSE ERROR on the runner (WP-43 r2
    # red: four suites failed to load). Guard the pattern itself.
    for f in sorted(pglob('dev/tests/test_*.m')):
        lines = open(f, encoding='utf-8', errors='replace').read().split('\n')
        for i, l in enumerate(lines):
            if l.strip().startswith('% tol:') and i > 0 \
                    and lines[i - 1].rstrip().endswith('...'):
                note('cat', '%s:%d: tol: comment inside a ... continuation '
                     '(parse error on load)' % (f, i + 1))

    # (d) tolerance provenance
    tol_re = re.compile(r"'(Abs|Rel)Tol'")
    debt = 0
    for f in sorted(pglob('dev/tests/test_*.m')):
        lines = open(f, encoding='utf-8', errors='replace').read().split('\n')
        for i, l in enumerate(lines):
            if tol_re.search(l):
                ctx = '\n'.join(lines[max(0, i - 3):i + 1])
                if 'tol:' not in ctx:
                    if f in TOL_SWEPT:
                        note('cat', '%s:%d: tolerance without a "tol:" '
                             'provenance comment (swept suite)' % (f, i + 1))
                    else:
                        debt += 1
    print('   unannotated tolerances: %d (must stay 0 -- full '
          'enforcement since WP-43 r2)' % debt)


def check_readme_sync():
    print('\n17. README agrees with the tree (doc sync gate, HANDOVER 6.2.4)')
    # WHY (WP-46, 2026-08-20): the README drifted for a month (v1-era test
    # counts and status) because no gate read it. Adopted principle (external
    # geoMap handover 2.7/2.8.2): one authority per fact, every other surface
    # CHECKED against it, never independently maintained.
    #   (a) the README title version must equal the Contents.m Version line
    #       (Contents.m is the ONE version source -- runMeta reads it too);
    #   (b) the README's suite/point-count claim ("N suite files, M test
    #       points") must equal the tree's counts (same static count as
    #       check 16).
    # L8 (WP-70): the counts claim lives in dev/DEVELOPMENT.md now; README
    # carries only intro / quick start / usage and keeps the version in its
    # title. Both stay CHECKED, neither independently maintained.
    readme = open('README.md', encoding='utf-8', errors='replace').read()
    devdoc = open('dev/DEVELOPMENT.md', encoding='utf-8',
                  errors='replace').read()
    contents = open('+lsqvce/Contents.m', encoding='utf-8',
                    errors='replace').read()
    mv = re.search(r'^%\s*Version\s+(\S+)', contents, re.M)
    mt = re.search(r'^#\s*lsqvce\s+v(\S+?)\s', readme, re.M)
    if not mv or not mt:
        note('readme', 'version line missing (Contents.m: %s, README title: '
             '%s)' % (bool(mv), bool(mt)))
    elif mv.group(1) != mt.group(1):
        note('readme', 'README title says v%s but Contents.m says %s'
             % (mt.group(1), mv.group(1)))

    files = sorted(pglob('dev/tests/test_*.m'))
    n_files = len(files)
    fn_re = re.compile(r'^function\s+(?:\w+\s*=\s*)?test\w+\s*\(')
    n_pts = 0
    for f in files:
        for l in open(f, encoding='utf-8', errors='replace'):
            if fn_re.match(l) and not re.match(r'function\s+tests\s*=', l):
                n_pts += 1
    mc = re.search(r'(\d+)\s+suite files,\s*\*{0,2}?(\d+)\s+test\s*\n?>?\s*points',
                   devdoc)
    if not mc:
        mc = re.search(r'(\d+)\s+suite files.{0,40}?(\d+)\s+test\s+points',
                       devdoc, re.S)
    if not mc:
        note('readme', 'dev/DEVELOPMENT.md carries no machine-checkable "N suite files, '
             'M test points" claim')
    else:
        cf, cp = int(mc.group(1)), int(mc.group(2))
        if cf != n_files or cp != n_pts:
            note('readme', 'dev/DEVELOPMENT.md claims %d files / %d points; the tree has '
                 '%d / %d' % (cf, cp, n_files, n_pts))
        else:
            print('   README claim %d files / %d points == tree; version v%s '
                  '== Contents.m' % (cf, cp, mv.group(1) if mv else '?'))


def check_tiers_and_case():
    print('\n18. Surface census: every public name has a tier row; no case-only '
          'file-name collisions (HANDOVER 11.9, WP-54)')
    # WHY (R4 WP-54, 2026-08-31): the public surface grew to 80 names with
    # no statement of who they are for. dev/docs_src/tiers.csv is the ONE tier
    # authority (USER / ENGINE / KERNEL); this check keeps it complete in
    # both directions so a new public name cannot land untiered and a
    # deleted one cannot linger. (b) guards the Windows working copy: two
    # names differing only by case (options.m / Options.m) cannot coexist on
    # a case-insensitive file system, so the tree must never contain them.
    import csv
    tiers_ok = {'USER', 'ENGINE', 'KERNEL'}
    path = 'dev/docs_src/tiers.csv'
    if not os.path.isfile(path):
        note('tiers', '%s missing' % path)
        return
    rows = list(csv.DictReader(open(path, encoding='utf-8')))
    names = {}
    for r in rows:
        if r.get('tier') not in tiers_ok:
            note('tiers', 'row %s has tier %r (expected one of %s)'
                 % (r.get('name'), r.get('tier'), sorted(tiers_ok)))
        if r['name'] in names:
            note('tiers', 'duplicate row for %s' % r['name'])
        names[r['name']] = r
    tree = set()
    for pkg, sub in (('', ''), ('engine.', '+engine/'), ('kernel.', '+kernel/'),
                     ('util.', '+util/')):
        for f in pglob('+lsqvce/%s*.m' % sub):
            n = os.path.basename(f)[:-2]
            if n != 'Contents':
                tree.add(pkg + n)
    for n in sorted(tree - set(names)):
        note('tiers', 'public name %s has no row in %s' % (n, path))
    for n in sorted(set(names) - tree):
        note('tiers', 'row %s names nothing in the tree' % n)
    by_tier = {}
    for r in rows:
        by_tier.setdefault(r['tier'], []).append(r['name'])
    print('   %d public names: USER %d, ENGINE %d, KERNEL %d' % (
        len(tree), len(by_tier.get('USER', [])),
        len(by_tier.get('ENGINE', [])), len(by_tier.get('KERNEL', []))))
    # (b) case-insensitive uniqueness, whole tree, .m and .py and .md
    seen = {}
    for root, _dirs, files in os.walk('.'):
        if '.git' in root or '__pycache__' in root:
            continue
        for fn in files:
            key = os.path.join(root, fn).lower()
            if key in seen:
                note('tiers', 'case-only file-name collision: %s vs %s '
                     '(breaks Windows working copies)'
                     % (os.path.join(root, fn), seen[key]))
            seen[key] = os.path.join(root, fn)




# ----------------------------------------------------------------------
# R6 WP-94 (2026-09-07): two doc-truth checks the docs level lacked. Found
# by measuring "are the docs in sync?": twelve help EXAMPLE blocks still
# called units deleted in 5.0 (adjust, adjustVCE(...).run()) -- check 13
# guards a fixed list of OLD names, so a deletion never reached it -- and
# thirteen live units had a page but no entry in html/index.html or
# helptoc.xml, so `doc` reached them by name only. Written with the
# assistance of Claude (Anthropic), model Claude Fable 5.1.

HELP_HEADING = re.compile(
    r'^(SYNTAX|INPUTS?|OUTPUTS?|OPTIONS?|NAME-VALUE.*|ERRORS?|EXAMPLES?.*|NOTES?|'
    r'REFERENCES?|ALGORITHM|CONVENTIONS?|WHY.*|WHAT.*|COST.*|SEE ALSO|See also.*|'
    r'BEST PRACTICES?|ACCURACY.*|LIMITATIONS?|CONTRACT.*|THE .*|PROPERTIES|METHODS|'
    r'WARNINGS?|IDENTIFIERS.*|RUNTIME.*|Created .*|Updated .*|Moved .*|Generated .*)')


def public_units():
    """{'fit', 'engine.makeLoader', ...} for every non-private, non-Contents
    unit file of the package (the names a help example may call)."""
    units = set()
    for f in pglob(PKG_GLOB, recursive=True):
        if '/private/' in f:
            continue
        u = _unit_of(f)
        if not u.endswith('Contents'):
            units.add(u)
    return units


def help_example_lines(path):
    """(lineno, text) of every comment line inside an EXAMPLE section of the
    file's help blocks (function help and method help alike)."""
    out = []
    in_ex = False
    for i, raw in enumerate(open(path, encoding='utf-8', errors='replace'), 1):
        s = raw.strip()
        if not s.startswith('%'):
            in_ex = False
            continue
        body = s.lstrip('%').strip()
        if HELP_HEADING.match(body):
            in_ex = body.upper().startswith('EXAMPLE')
            continue
        if in_ex:
            out.append((i, body))
    return out


def check_help_examples_call_live_units():
    print('\n24. help EXAMPLE blocks call units that exist')
    units = public_units()
    n = 0
    for f in sorted(pglob(PKG_GLOB, recursive=True)):
        if '/private/' in f:
            continue
        for i, text in help_example_lines(f):
            for u in re.findall(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Za-z]\w*)', text):
                n += 1
                if u not in units:
                    note('example', '%s:%d example calls lsqvce.%s, which does not exist'
                         % (f, i, u))
    print('   %d unit references in EXAMPLE blocks checked against %d units'
          % (n, len(units)))


def check_index_lists_every_unit():
    print('\n25. every unit page is reachable from html/index.html and helptoc.xml')
    index = open('html/index.html', encoding='utf-8', errors='replace').read()
    toc = open('helptoc.xml', encoding='utf-8', errors='replace').read()
    units = sorted(public_units())
    for u in units:
        page = 'lsqvce.%s.html' % u
        if page not in index:
            note('index', 'html/index.html has no row for lsqvce.%s' % u)
        if page not in toc:
            note('index', 'helptoc.xml has no entry for lsqvce.%s' % u)
    print('   %d units checked' % len(units))
MD_SCOPE = ['README.md', 'CLAUDE.md', 'RELEASE_NOTES.md', 'dev/DEVELOPMENT.md',
            'dev/SETUP_GUIDE.md', 'dev/GRAVFIRE_INTERFACE.md',
            'dev/references/README.md']
MD_EXT = ('.m', '.md', '.html', '.pdf', '.py', '.csv', '.xml', '.json', '.sh',
          '.yml', '.mat')
MD_SKIP = ('http', 'mailto', '#', 'dev/references/local/')
MD_WILD = ('*', '?', '<', '>', '...', '$', '{')
_MD_LINK = re.compile(r'\[[^\]]*\]\(([^)]+)\)')
_MD_TICK = re.compile(r'`([^`\n]+)`')


def _md_pointers(text):
    """Every repository pointer in a markdown text: (target, kind)."""
    for m in _MD_LINK.finditer(text):
        t = m.group(1).split('#')[0].strip()
        if t:
            yield t, 'link'
    for m in _MD_TICK.finditer(text):
        raw = m.group(1).strip()
        if not raw:
            continue
        # a backticked span is often a COMMAND ("python3 dev/docs_src/x.py"):
        # the pointer is its last token, and trailing prose punctuation is not
        # part of a path
        t = raw.split()[-1].rstrip('.,;:')
        if '/' in t and t.endswith(MD_EXT):
            yield t, 'tick'


def check_md_pointers():
    print('\n26. markdown pointers resolve (R6 WP-99)')
    # WHY: check 5 reads links inside html/ only, so NOTHING looked at a
    # pointer in the repository's own prose. MEASURED when this was written:
    # dev/DEVELOPMENT.md still sent the reader to four HTML manuals deleted in
    # the docs level of 2026-09-04, plus a build script and a benchmark that no
    # longer exist -- ten dead pointers that checks 5, 24 and 25 all walked
    # past, in the document a developer reads FIRST.
    #
    # SCOPE is deliberate, not everything with a .md suffix. These are the
    # documents that TELL A READER WHERE TO GO. The records -- dev/HANDOVER.md
    # and dev/HISTORY_v1.md -- are excluded on purpose: their job is to say what
    # WAS true, and half their value is naming files that were deleted.
    # dev/docs_src/manual/*.md is excluded because check 23 already rebuilds
    # the Manual and compares it byte for byte.
    n = 0
    for f in MD_SCOPE:
        if not os.path.exists(f):
            continue
        text = open(f, encoding='utf-8', errors='replace').read()
        for t, kind in _md_pointers(text):
            if t.startswith(MD_SKIP) or any(w in t for w in MD_WILD):
                continue
            n += 1
            near = os.path.normpath(os.path.join(os.path.dirname(f), t))
            if os.path.exists(t) or os.path.exists(near):
                continue
            note('md-pointer', '%s points at %s (%s), which does not exist'
                 % (f, t, kind))
    print('   %d pointers checked in %d documents'
          % (n, sum(1 for f in MD_SCOPE if os.path.exists(f))))


def check_options_are_consumed():
    print('\n27. every public Options property has a consumer (R6 WP-101)')
    # WHY: MinPasses and LogFile shipped for a whole review as DEAD KNOBS --
    # each appeared exactly three times in the package (the help line, the
    # property declaration, the disp() group table) and nowhere else. A
    # documented option that nothing reads is worse than a missing one: the
    # user sets it, nothing refuses it, and the run quietly ignores what they
    # asked for. No other check could see it -- the help was well-formed, the
    # page matched the help, the identifier catalogue was closed.
    #
    # THE RULE: a property is CONSUMED if its name appears either inside
    # Options.engineArgs (the ONE hand-off to the drivers, which may forward it
    # under a different engine name) or anywhere in the package outside
    # Options.m. Anything else is dead surface.
    src = open('+lsqvce/Options.m', encoding='utf-8', errors='replace').read()
    mprop = re.search(r'^    properties\s*$(.*?)^    end\s*$', src, re.S | re.M)
    if not mprop:
        note('options', '+lsqvce/Options.m: no properties block found')
        return
    props = re.findall(r'^        ([A-Z]\w*)\s', mprop.group(1), re.M)
    marg = re.search(r'function args = engineArgs\(.*?\n(.*?)\n        end', src, re.S)
    handoff = marg.group(1) if marg else ''
    elsewhere = ''
    for f in pglob(PKG_GLOB, recursive=True):
        if f.replace('\\', '/').endswith('+lsqvce/Options.m'):
            continue
        elsewhere += open(f, encoding='utf-8', errors='replace').read()
    for p in props:
        if re.search(r'\b%s\b' % re.escape(p), handoff):
            continue
        if re.search(r'\b%s\b' % re.escape(p), elsewhere):
            continue
        note('options', 'Options.%s is declared and documented but NOTHING reads '
             'it: not engineArgs, not anywhere else in +lsqvce' % p)
    print('   %d public properties checked' % len(props))


def _help_descriptive(path):
    """The help block of an .m file up to (excluding) its 'See also' line --
    the part that makes claims; the stamps after it are history."""
    lines = open(path, encoding='utf-8', errors='replace').read().split('\n')
    out = []
    started = False
    for ln in lines:
        st = ln.strip()
        if not started:
            if st.startswith('%'):
                started = True
            elif st and not st.startswith('function') and not st.startswith('classdef'):
                break
            else:
                continue
        if started:
            if not st.startswith('%'):
                break
            if re.match(r'%\s*See also', st, re.I):
                break
            out.append(ln)
    return '\n'.join(out)


def check_superseded_claims():
    print('\n28. no live surface still makes a superseded claim (R6 WP-106)')
    # WHY: check 22 proved every page equals its help block and check 23 that
    # the Manual equals its sources -- and a sentence that was TRUE at 5.2 and
    # FALSE at 6.1 passed both, twice in one day: Options.m still offered
    # Noise.Grid "auto" after WP-101 retired it, and Manual chapter 4 still
    # said Robust and Estimate="noise" were "not combinable" after WP-105 made
    # them run. Byte-equality cannot see a claim go stale. This check can see
    # ONE thing: a phrase that once described a behaviour that no longer
    # exists. Every retirement adds its phrase to superseded_claims.csv, and
    # from then on it cannot come back anywhere it would mislead a reader.
    csv_path = os.path.join('dev', 'docs_src', 'superseded_claims.csv')
    rows = []
    for ln in open(csv_path, encoding='utf-8'):
        ln = ln.rstrip('\n')
        if not ln.strip() or ln.lstrip().startswith('#'):
            continue
        parts = ln.split(',', 2)
        if len(parts) < 3:
            note('superseded', 'superseded_claims.csv: malformed row %r' % ln)
            continue
        rows.append((parts[0].strip(), parts[1].strip(), parts[2].strip()))
    surfaces = []
    for f in pglob(PKG_GLOB, recursive=True):
        surfaces.append((f, _help_descriptive(f)))
    for pat in ('dev/docs_src/manual/*.md', 'README.md', 'dev/GRAVFIRE_INTERFACE.md',
                'dev/DEVELOPMENT.md', 'template/*.m', 'examples/*.m'):
        for f in pglob(pat):
            surfaces.append((f, open(f, encoding='utf-8', errors='replace').read()))
    n = 0
    for pattern, by, claim in rows:
        rx = re.compile(pattern)
        for f, text in surfaces:
            for i, ln in enumerate(text.split('\n'), 1):
                if rx.search(ln):
                    note('superseded', '%s:%d still makes a claim retired by %s (%s): %r'
                         % (f, i, by, claim, ln.strip()[:100]))
                    n += 1
    print('   %d superseded claims checked on %d surfaces; hits: %d' % (len(rows), len(surfaces), n))


def main():
    if not os.path.isdir('+lsqvce'):
        sys.exit('run from the toolbox root (the folder holding +lsqvce)')
    print('lsqvce tree audit')
    print('=' * 60)
    check_test_registration()
    check_removed_names()
    check_signatures()
    check_doc_coverage()
    check_links()
    check_no_js()
    check_balance()
    check_banners()
    check_spec_fields()
    check_duplicate_functions()
    check_dedup_lock()
    check_deprecated_words()
    check_addpath_targets()
    check_html_examples()
    check_error_ids()
    check_h1_names()
    check_test_categories()
    check_pages_match_help()
    check_manual_built_and_complete()
    check_readme_sync()
    check_tiers_and_case()
    check_catch_ratchet()
    check_id_catalogue()
    check_layering_ratchet()
    check_help_examples_call_live_units()
    check_index_lists_every_unit()
    check_md_pointers()
    check_options_are_consumed()
    check_superseded_claims()
    check_manual_pdf_version()
    print('\n' + '=' * 60)
    if FINDINGS:
        print('%d finding(s).' % len(FINDINGS))
        return 1
    print('CLEAN - no findings.')
    return 0




# ----------------------------------------------------------------------
# R5 Level M0 (WP-78, 2026-09-01): three RATCHETS and a catalogue. A ratchet
# lists the debt that exists today in a baseline file, refuses any NEW
# instance and refuses a STALE row, so the debt can only shrink and every
# shrink is a committed edit of the baseline. Written with the assistance of
# Claude (Anthropic), model Claude Fable 5.1.

PKG_GLOB = '+lsqvce/**/*.m'


def _unit_of(path):
    """'+lsqvce/+engine/foo.m' -> 'engine.foo'; root files -> 'foo'."""
    u = path.replace('+lsqvce/', '').replace('.m', '')
    u = u.replace('+util/', 'util.').replace('+engine/', 'engine.')
    u = u.replace('+kernel/', 'kernel.').replace('private/', 'private.')
    return u


def _read_baseline(path, ncol):
    rows = []
    if not os.path.isfile(path):
        return rows
    for line in open(path, encoding='utf-8'):
        s = line.strip()
        if not s or s.startswith('#') or s.startswith('site,') or s.startswith('id,') \
                or s.startswith('caller,'):
            continue
        parts = [p.strip() for p in s.split(',', ncol - 1)]
        if len(parts) == ncol:
            rows.append(tuple(parts))
    return rows


def bare_catch_sites():
    """Return [(site, kind)] for every catch block that is EMPTY (no code at
    all) or BARE (no identifier bound AND nothing inspected: `catch` with a
    body that never rethrows, records or names an id). Site = unit:line."""
    sites = []
    for f in sorted(pglob(PKG_GLOB, recursive=True)):
        lines = open(f, encoding='utf-8', errors='replace').read().split('\n')
        for i, raw in enumerate(lines):
            m = re.match(r'\s*catch\b\s*(\w+)?\s*(%.*)?$', raw)
            if not m:
                continue
            bound = m.group(1)
            body = []
            depth = 0
            for j in range(i + 1, len(lines)):
                s = strip_matlab(lines[j])
                if re.match(r'\s*end\b', s) and depth == 0:
                    break
                if re.match(r'\s*(if|for|while|try|switch|parfor)\b', s):
                    depth += 1
                elif re.match(r'\s*end\b', s):
                    depth -= 1
                if s.strip():
                    body.append(s.strip())
            if not body:
                sites.append(('%s:%d' % (_unit_of(f), i + 1), 'empty'))
            elif not bound:
                sites.append(('%s:%d' % (_unit_of(f), i + 1), 'bare'))
    return sites


def check_catch_ratchet():
    print('\n19. catch-block ratchet (R5 M0): no new empty/bare catch')
    # WHY: five empty and thirteen bare `catch` blocks (R5 finding H-07)
    # convert an undefined variable or a shape error into a silently
    # different run -- the fourth outcome the error doctrine forbids
    # (theory_r5 ch. 7). M1 rewrites every site; this check stops the count
    # from growing in the meantime and refuses a baseline row that no longer
    # exists, so the file is always the truth.
    base = dict(_read_baseline('dev/docs_src/catch_baseline.csv', 2))
    now = dict(bare_catch_sites())
    print('   empty/bare catch sites: %d (baseline %d; doctrine target 0)'
          % (len(now), len(base)))
    for site, kind in sorted(now.items()):
        if site not in base:
            note('catch', 'NEW %s catch at %s -- refuse/report/proceed; '
                          'never a bare default (theory_r5 ch. 7)' % (kind, site))
    for site in sorted(base):
        if site not in now:
            note('catch', 'stale baseline row %s (fixed or moved): delete '
                          'it from dev/docs_src/catch_baseline.csv' % site)


RAISE_PATS = (r"error\(\s*'([\w]+:[\w:]+)'",
              r"MException\(\s*'([\w]+:[\w:]+)'",
              r"warning\(\s*'([\w]+:[\w:]+)'",
              r"log\.warn\w*\(\s*'([\w]+:[\w:]+)'",
              r"warnOnce\(\s*'([\w]+:[\w:]+)'",
              r"'([\w]+:[\w:]+)'\s*,\s*\.\.\.")


def raised_ids():
    """{id: set(units)} for every lsqvce:* identifier the PACKAGE raises."""
    out = {}
    for f in sorted(pglob(PKG_GLOB, recursive=True)):
        src = open(f, encoding='utf-8', errors='replace').read()
        src = src.replace('("', "('").replace('",', "',")
        code = '\n'.join(s for _, s in code_lines(f))
        # raise sites live in code; the string-stripped text loses the id,
        # so scan the raw source but only lines that survive as code
        raw_lines = src.split('\n')
        codeset = {k for k, _ in code_lines(f)}
        text = '\n'.join(l for k, l in enumerate(raw_lines, 1) if k in codeset)
        for pat in RAISE_PATS:
            for i in re.findall(pat, text):
                if i.startswith('lsqvce:'):
                    out.setdefault(i, set()).add(_unit_of(f))
        del code
    return out


def check_id_catalogue():
    print('\n20. error-identifier catalogue (R5 M0): closed both ways')
    # WHY: 369 ids in 40+ unit namespaces with no authority (finding H-09).
    # The catalogue dev/docs_src/error_ids.csv IS the authority from M0 on:
    # an id raised but not catalogued is a finding; an id catalogued but no
    # longer raised is a finding (so deletions in M2 are recorded, not
    # forgotten). Documented-in-help is the second obligation: reported as
    # the debt number here, enforced in M1.
    cat = {}
    for row in _read_baseline('dev/docs_src/error_ids.csv', 4):
        cat[row[0]] = row[1:]
    raised = raised_ids()
    print('   raised lsqvce:* identifiers: %d; catalogued: %d'
          % (len(raised), len(cat)))
    for i in sorted(raised):
        if i not in cat:
            note('id-catalogue', '%s is raised by %s but not in '
                 'dev/docs_src/error_ids.csv' % (i, ','.join(sorted(raised[i]))))
    for i in sorted(cat):
        if i not in raised:
            note('id-catalogue', '%s is catalogued but no package code '
                 'raises it (delete the row or restore the raise)' % i)
    # documented in the raising unit's own help text?
    undocumented = 0
    help_text = {}
    for f in sorted(pglob(PKG_GLOB, recursive=True)):
        help_text[_unit_of(f)] = '\n'.join(
            l for l in open(f, encoding='utf-8', errors='replace')
            if l.lstrip().startswith('%'))
    for i, units in raised.items():
        if not any(i in help_text.get(u, '') for u in units):
            undocumented += 1
    print('   raised ids not documented in their unit\'s help: %d '
          '(M1 debt; the doctrine target is 0)' % undocumented)


TIER_RANK = {'USER': 3, 'ENGINE': 2, 'KERNEL': 1, 'UTIL': 0}


def layering_edges():
    """Return {(caller, callee, callerTier, calleeTier)} for every UPWARD
    qualified call in the package (callee tier higher than caller tier).
    Tiers from dev/docs_src/tiers.csv; util.* is its own tier UTIL."""
    tiers = {}
    if os.path.isfile('dev/docs_src/tiers.csv'):
        for row in open('dev/docs_src/tiers.csv', encoding='utf-8'):
            parts = row.strip().split(',', 2)
            if len(parts) >= 2 and parts[0] != 'name':
                name, tier = parts[0].strip(), parts[1].strip()
                if name.startswith('util.'):
                    tier = 'UTIL'
                tiers[name] = tier
    edges = set()
    for f in sorted(pglob(PKG_GLOB, recursive=True)):
        if '/private/' in f:
            continue
        caller = _unit_of(f)
        ct = tiers.get(caller)
        if ct is None:
            continue
        text = '\n'.join(s for _, s in code_lines(f))
        for m in re.finditer(r'lsqvce\.((?:util|engine|kernel)\.\w+|[A-Z]\w+|fit|report)', text):
            callee = m.group(1)
            if callee == caller:
                continue
            kt = tiers.get(callee)
            if kt is None:
                continue
            if TIER_RANK[kt] > TIER_RANK[ct]:
                edges.add((caller, callee, ct, kt))
    return edges


def check_layering_ratchet():
    print('\n21. layering ratchet (R5 M0): calls go downward only')
    # WHY (CODING_GUIDE Part 3): a layering rule nothing checks is a
    # diagram. The engine constructs lsqvce.Results and calls lsqvce.report
    # today (USER-tier names) -- known, listed, and the M4/M5 redesign
    # removes them. Until then: no NEW upward edge, no stale baseline row.
    base = {(r[0], r[1]) for r in _read_baseline('dev/docs_src/layering_baseline.csv', 4)}
    now = layering_edges()
    print('   upward edges: %d (baseline %d; target 0 after M5)'
          % (len(now), len(base)))
    nowk = {(e[0], e[1]) for e in now}
    for e in sorted(now):
        if (e[0], e[1]) not in base:
            note('layering', 'NEW upward call %s (%s) -> %s (%s)'
                 % (e[0], e[2], e[1], e[3]))
    for c, k in sorted(base):
        if (c, k) not in nowk:
            note('layering', 'stale baseline row %s -> %s: delete it from '
                 'dev/docs_src/layering_baseline.csv' % (c, k))


def write_baselines():
    """Regenerate the three baseline/catalogue files from the tree (M0
    bootstrap, and the documented way to accept a deliberate change:
    `python3 dev/docs_src/audit_tree.py --write-baselines`, then READ the
    diff -- a ratchet that is regenerated blindly is not a ratchet)."""
    stamp = 'generated by audit_tree.py --write-baselines'
    with open('dev/docs_src/catch_baseline.csv', 'w', encoding='utf-8') as fh:
        fh.write('# Empty/bare catch sites accepted as DEBT (R5 M0, check 19).\n'
                 '# One row per site (unit:line, kind). M1 (WP-80) drives this to\n'
                 '# zero; a fixed site is DELETED here in the same commit.\n'
                 '# %s\nsite,kind\n' % stamp)
        for site, kind in sorted(bare_catch_sites()):
            fh.write('%s,%s\n' % (site, kind))
    with open('dev/docs_src/error_ids.csv', 'w', encoding='utf-8') as fh:
        fh.write('# Error-identifier catalogue (R5 M0, check 20) -- THE authority\n'
                 '# for lsqvce:* identifiers: closed both ways against the tree.\n'
                 '# id, units raising it (;-separated), kind (error|warning|mixed), since\n'
                 '# %s\nid,units,kind,since\n' % stamp)
        raised = raised_ids()
        kinds = {}
        for f in sorted(pglob(PKG_GLOB, recursive=True)):
            src = open(f, encoding='utf-8', errors='replace').read()
            src = src.replace('("', "('").replace('",', "',")
            for pat in RAISE_PATS[:2]:
                for i in re.findall(pat, src):
                    kinds.setdefault(i, set()).add('error')
            for pat in RAISE_PATS[2:5]:
                for i in re.findall(pat, src):
                    kinds.setdefault(i, set()).add('warning')
        for i in sorted(raised):
            k = kinds.get(i, {'error'})
            kind = 'mixed' if len(k) > 1 else next(iter(k))
            fh.write('%s,%s,%s,v4.0\n' % (i, ';'.join(sorted(raised[i])), kind))
    with open('dev/docs_src/layering_baseline.csv', 'w', encoding='utf-8') as fh:
        fh.write('# Upward calls accepted as DEBT (R5 M0, check 21): caller, callee,\n'
                 '# callerTier, calleeTier. M4/M5 remove them; a removed edge is\n'
                 '# DELETED here in the same commit.\n# %s\ncaller,callee,callerTier,calleeTier\n' % stamp)
        for e in sorted(layering_edges()):
            fh.write('%s,%s,%s,%s\n' % e)
    print('baselines written (read the diff before committing)')


def check_pages_match_help():
    print('\n22. every unit page equals a regeneration from its help block (no doc drift)')
    # R5 docs level (2026-09-04): the API pages are DERIVED from the help blocks
    # (build_missing_pages.render); a page that differs from a fresh render is
    # drift. Regenerate with `python3 dev/docs_src/build_missing_pages.py all`.
    import importlib.util as _il
    spec = _il.spec_from_file_location('bmp', os.path.join('dev', 'docs_src', 'build_missing_pages.py'))
    bmp = _il.module_from_spec(spec)
    spec.loader.exec_module(bmp)
    drifted = bmp.main('check')
    for name in drifted:
        note('docdrift', 'html/lsqvce.%s.html differs from a regeneration of its help block '
             '(run build_missing_pages.py all)' % name)
    print('   pages checked; drifted: %d' % len(drifted))


def check_manual_built_and_complete():
    print('\n23. the Manual equals a fresh build and names every Options knob and group field')
    # R5 docs level (2026-09-04): html/lsqvce_manual.html is built by
    # dev/docs_src/build_manual.py from dev/docs_src/manual/*.md plus generated
    # chapters (API from the help blocks, tiers from Contents.m, migration from
    # RELEASE_NOTES.md, identifiers from error_ids.csv). Stale = drift.
    import importlib.util as _il
    spec = _il.spec_from_file_location('bm', os.path.join('dev', 'docs_src', 'build_manual.py'))
    bm = _il.module_from_spec(spec)
    try:
        spec.loader.exec_module(bm)
        if bm.main(check=True):
            note('manual', 'html/lsqvce_manual.html differs from a fresh build (run build_manual.py)')
    except Exception as e:  # pandoc missing etc.: report, do not hide
        note('manual', 'the Manual could not be rebuilt for comparison: %s' % e)
    src = open(os.path.join('dev', 'docs_src', 'manual', '04_options.md'), encoding='utf-8').read()
    opt = open(os.path.join('+lsqvce', 'Options.m'), encoding='utf-8').read()
    i = opt.index('    properties'); j = opt.index('    end', i)
    props = re.findall(r'^\s{8}([A-Z]\w*)\s', opt[i:j], re.M)
    tail = opt[opt.index('methods (Static)'):opt.index('methods (Static, Access = private)')]
    fields = set(re.findall(r"'([A-Z][A-Za-z0-9]*)', ", tail))
    for p in props:
        if not re.search(r'`%s`' % p, src):
            note('manual', 'Options.%s is not described in the Manual chapter 4' % p)
    for f in sorted(fields):
        if f in props: continue
        if not re.search(r'`%s`' % f, src):
            note('manual', 'the Options group field %s is not described in the Manual chapter 4' % f)
    print('   %d knobs and %d group fields checked' % (len(props), len(fields - set(props))))


# ----------------------------------------------------------------------
# 29 (R6 WP-109, 2026-09-20): the PDF Manual's version surface. Written with
# the assistance of Claude (Anthropic), model Claude Opus 5.

def pdf_doc_title(path):
    """The DOCUMENT title of a PDF (the /Title of its Info dictionary), or None.

    Why this is hand-rolled rather than shelled out to pdftotext/pypdf: the
    audit must run on the CI runner with nothing but Python, and it needs ONE
    string, not a text layer. A pandoc/xelatex PDF hides its Info dictionary
    inside a FlateDecode'd object stream, so the bytes are searched after
    inflating every stream. The dictionary is identified by /Producer (every
    TeX PDF in this tree has one) and delimited by its own << >>, so a
    BOOKMARK /Title -- there are 51 of them in the Manual -- cannot be
    mistaken for the document's.

    Both PDF string forms are handled, because both occur in this tree:
    lsqvce_manual.pdf (xelatex/xdvipdfmx) writes UTF-16BE hex <feff...>, and
    lsqvce_theory.pdf (pdfTeX/hyperref) writes an octal-escaped literal.
    """
    raw = open(path, 'rb').read()
    bufs = [raw]
    for m in re.finditer(rb'stream\r?\n', raw):
        a = m.end()
        b = raw.find(b'endstream', a)
        if b < 0:
            continue
        try:
            bufs.append(zlib.decompress(raw[a:b]))
        except Exception:
            pass                       # not a Flate stream: nothing to read
    titles = []
    for buf in bufs:
        for m in re.finditer(rb'/Producer', buf):
            lo = buf.rfind(b'<<', 0, m.start())
            hi = buf.find(b'>>', m.start())
            if lo < 0 or hi < 0:
                continue
            info = buf[lo:hi]
            h = re.search(rb'/Title\s*<(feff[0-9a-fA-F]+)>', info)
            if h:
                titles.append(binascii.unhexlify(h.group(1))
                              .decode('utf-16-be').lstrip('\ufeff'))
                continue
            t = re.search(rb'/Title\s*\((.*?)(?<!\\)\)', info, re.S)
            if t:
                titles.append(pdf_literal(t.group(1)))
    uniq = sorted(set(titles))
    return uniq[0] if len(uniq) == 1 else None


def pdf_literal(s):
    """Decode a PDF literal string: octal escapes, then UTF-16BE if it has a BOM."""
    out = bytearray()
    i = 0
    while i < len(s):
        if s[i:i + 1] == b'\\' and s[i + 1:i + 2].isdigit():
            out.append(int(s[i + 1:i + 4], 8))
            i += 4
        elif s[i:i + 1] == b'\\':
            out.append(s[i + 1])
            i += 2
        else:
            out.append(s[i])
            i += 1
    b = bytes(out)
    if b[:2] == b'\xfe\xff':
        return b[2:].decode('utf-16-be', 'replace')
    return b.decode('latin-1')


def check_manual_pdf_version():
    print('\n29. the PDF Manual carries the Contents.m version in its title')
    # WHY (R6 WP-109, out of WP-108): lsqvce_manual.pdf sat at 6.1.0 while the
    # tree was at 6.1.1 and NOTHING read it. Check 23 rebuilds the HTML Manual
    # and compares it byte for byte; the PDF is built from the same markdown by
    # the same script but cannot be compared that way -- a rebuild of the SAME
    # source differs in the font subset tags, the document /ID and the xref
    # offsets, which xdvipdfmx randomises per run (measured: 21 of 141 streams,
    # every difference accounted for, none of them content). So this check
    # reads the ONE fact that must not drift: the title build_manual.py writes
    # from Contents.m, the same single version source check 17 (a) holds the
    # README title to.
    #
    # lsqvce_theory.pdf is deliberately NOT checked: its title is
    # "lsqvce - Theory, practice, ..." and carries no version at all. It is a
    # dated record and WP-106 gave it banners instead; there is no version
    # surface here to gate.
    #
    # It FAILS CLOSED. A missing file, an unreadable Info dictionary or two
    # different answers are findings, not silence -- the instrument saying
    # nothing must never be read as the artefact being right.
    pdf = 'lsqvce_manual.pdf'
    if not os.path.isfile(pdf):
        note('manual-pdf', '%s is missing (build_manual.py writes it)' % pdf)
        return
    contents = open(os.path.join('+lsqvce', 'Contents.m'),
                    encoding='utf-8', errors='replace').read()
    mv = re.search(r'^%\s*Version\s+(\S+)', contents, re.M)
    if not mv:
        note('manual-pdf', '+lsqvce/Contents.m has no Version line to check against')
        return
    title = pdf_doc_title(pdf)
    if title is None:
        note('manual-pdf', '%s: no single document title could be read from the '
             'Info dictionary (the check cannot pass on an unreadable PDF)' % pdf)
        return
    want = 'lsqvce %s \u2014 User Manual' % mv.group(1)
    if title != want:
        note('manual-pdf', '%s title is "%s"; Contents.m says %s, so it must be '
             '"%s" (run build_manual.py)' % (pdf, title, mv.group(1), want))
    else:
        print('   %s == Contents.m %s' % (title, mv.group(1)))


if __name__ == '__main__':
    if '--write-baselines' in sys.argv:
        if not os.path.isdir('+lsqvce'):
            sys.exit('run from the toolbox root (the folder holding +lsqvce)')
        write_baselines()
        sys.exit(0)
    sys.exit(main())
