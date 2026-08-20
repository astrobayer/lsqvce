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

USAGE
    python3 docs_src/audit_tree.py            # from the toolbox root
    echo $?                                   # 0 = clean, 1 = findings

Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""

import glob
import os
import re
import sys
import xml.etree.ElementTree as ET

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
    runner = os.path.join('tests', 'run_all_lsqvce.m')
    if not os.path.exists(runner):
        note('tests', 'tests/run_all_lsqvce.m is missing')
        return
    reg = set(re.findall(r"'(test_\w+)\.m'", open(runner, encoding='utf-8').read()))
    disk = {os.path.basename(f)[:-2]
            for f in glob.glob('tests/*.m') + glob.glob('tests/accum/*.m')
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
        for f in glob.glob('tests/**/*.m', recursive=True):
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
    contents = ''
    if os.path.exists('+lsqvce/Contents.m'):
        contents = open('+lsqvce/Contents.m', encoding='utf-8').read()
    funcs = []
    for f in sorted(glob.glob('+lsqvce/*.m')):
        b = os.path.basename(f)[:-2]
        if b != 'Contents':
            funcs.append((b, 'html/lsqvce.%s.html' % b))
    for f in sorted(glob.glob('+lsqvce/+util/*.m')):
        b = os.path.basename(f)[:-2]
        funcs.append(('util.' + b, 'html/lsqvce.util.%s.html' % b))
    for name, page in funcs:
        if not os.path.exists(page):
            note('docs', 'no html page for %s (doc lsqvce.%s fails)' % (name, name))
        short = name.split('.')[-1]
        if not re.search(r'\b%s\b' % re.escape(short), contents):
            note('docs', '%s is not listed in Contents.m' % name)
    print('   %d functions checked' % len(funcs))

    for page in sorted(glob.glob('html/lsqvce.*.html')):
        n = os.path.basename(page)[len('lsqvce.'):-len('.html')]
        m = ('+lsqvce/+util/%s.m' % n[len('util.'):]) if n.startswith('util.') \
            else ('+lsqvce/%s.m' % n)
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
    for f in glob.glob('html/*.html'):
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
    for f in glob.glob('html/**/*.html', recursive=True):
        s = open(f, encoding='utf-8', errors='replace').read()
        if '<script' in s.lower():
            note('js', '%s contains a <script> tag' % f)
            n += 1
    for f in glob.glob('**/*.js', recursive=True):
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
    for f in glob.glob('+lsqvce/**/*.m', recursive=True):
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
            raw = open(f, encoding='utf-8', errors='replace').read()
            code = '\n'.join(l for l in raw.split('\n')
                             if not l.lstrip().startswith('%'))
            for m in pat.finditer(code):
                hits.setdefault(m.group(1), []).append(f)
    ok = True
    for name in ('makeKeyFcn', 'arcBatchWhiten'):
        where = hits.get(name, [])
        want = os.path.join('+lsqvce', 'private', name + '.m')
        if len(where) != 1 or os.path.normpath(where[0]) != os.path.normpath(want):
            note('dedup-lock', '%s must be defined EXACTLY once, in %s (found: %s)'
                 % (name, want, ', '.join(where) or 'nowhere'))
            ok = False
    if ok:
        print('    OK - both defined once, under +lsqvce/private/')

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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
        for f in glob.glob(d + '/**/*.m', recursive=True):
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
    for f in glob.glob('*.m') + glob.glob('tests/*.m'):
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
        'vceGroupsBy':  'groupComponentsBy',
        r'\.shape\(':   '.stochasticModel(',
        r'\.scale\(':   '.varianceComponents(',
        r'adjustVCE\(nc,':    'adjustVCE(systems, days) -- nc is derived',
        r'adjustBlock\(nc,':  'adjustBlock(systems, days) -- nc is derived',
        r'adjustBlock\(\d+,': 'adjustBlock(systems, days) -- nc is derived',
        r'P\.convergence\(':  'R.convergence() -- moved to Results',
        r'runtests\(':        'run_all_lsqvce',
    }
    # Only <pre> blocks: prose may legitimately DISCUSS an old name (e.g. a
    # migration note), but an example is what gets copied.
    pre = re.compile(r'<pre\b[^>]*>(.*?)</pre>', re.S)
    for f in sorted(glob.glob('html/*.html')):
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
    for f in glob.glob('+lsqvce/**/*.m', recursive=True) + \
             glob.glob('template/*.m') + glob.glob('tests/*.m'):
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
    for f in glob.glob('+lsqvce/**/*.m', recursive=True) + \
             glob.glob('html/*.html') + glob.glob('template/*.m'):
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
        for f in sorted(glob.glob(d + '/**/*.m', recursive=True)):
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
    cats_ok = {'contract', 'reference', 'precision', 'speed',
               'robustness', 'vectorisation'}
    tags_ok = {'diagnostics', 'provenance'}
    TOL_SWEPT = ['tests/test_signalCore.m', 'tests/test_precision.m']

    decl_re = re.compile(
        r'^\s*% categories:\s*([a-z, ]+?)\s*(?:\|\s*tags:\s*([a-z, ]+?)\s*)?$')
    fn_re = re.compile(r'function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(')

    unit_cat = {}          # unit -> set of categories reaching it
    n_points = 0
    for f in sorted(glob.glob('tests/test_*.m')):
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
            for u in set(re.findall(r'lsqvce\.(util\.\w+|[A-Za-z]\w+)', body)):
                unit_cat.setdefault(u, set()).update(cats)

    units = []
    for f in sorted(glob.glob('+lsqvce/**/*.m', recursive=True)):
        if '/private/' in f:
            continue
        u = f.replace('+lsqvce/', '').replace('.m', '').replace('+util/',
                                                                'util.')
        if u != 'Contents':
            units.append(u)

    exempt = {}
    if os.path.isfile('tests/category_exemptions.csv'):
        for row in open('tests/category_exemptions.csv',
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
    print('   category gaps in the matrix (informational): %d' % gaps)

    # (d) tolerance provenance
    tol_re = re.compile(r"'(Abs|Rel)Tol'")
    debt = 0
    for f in sorted(glob.glob('tests/test_*.m')):
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
    print('   unannotated tolerances outside the swept suites: %d '
          '(debt, drive to 0 suite by suite)' % debt)


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
    print('\n' + '=' * 60)
    if FINDINGS:
        print('%d finding(s).' % len(FINDINGS))
        return 1
    print('CLEAN - no findings.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
