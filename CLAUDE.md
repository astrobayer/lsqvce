# CLAUDE.md — pointer, not a status document

`dev/HANDOVER.md` is the ONLY status-and-plans document
(SSOT). Binding for every session: **§0** (standing rules + the loop),
**§1** (the work-package ledger — the only place a status lives),
**§6** (test/doc guidelines), and **§12** (the current review and its round
notes — §12.14 opened R6, the last subsection is the round in progress).
§8–§11 are the closed reviews R2…R4: history, not instructions.
`dev/HISTORY_v1.md` is frozen v1 history. Do not create a third document.

## Where the project stands

Read it, do not restate it: the version is `+lsqvce/Contents.m` (the ONE
version source, which the README title and `runMeta` are checked against),
the status of every work package is §1, and what the round in progress
decided is the last §12 subsection. Nothing about status is recorded in this
file — a second copy would be a second authority.

## The three commands

```matlab
lsqvce_setup                    % USER path (root + examples) + Help docs
cd('dev/tests'); ci_gate        % THE gate: full suite + log, or error
                                % (manages its path itself, restores it)
```

```bash
python3 dev/docs_src/audit_tree.py     # static audit, exit 0 required
python3 dev/docs_src/audit_selftest.py # the audit's own fault injection
```

The audit runs on any platform since WP-93 (`pglob`), with ONE exception:
check 23 rebuilds the Manual and compares it byte for byte against the
committed build, so it needs the same pandoc as the CI runner. On a Windows
box without pandoc the audit's honest baseline is **exactly one finding**
("the Manual could not be rebuilt"); anything more is yours. To rebuild the
Manual locally, use the Debian pandoc under WSL — CI is the authority.

The **postflight** cited in every PR is `ci_gate` green — and that is now the
whole of it. Since R6 WP-99 the examples are IN the gate: `test_examples.m`
runs every script in `examples/` once, in isolation, and refuses to let one be
added without a row in its call table. Before that the example run lived in a
scratch directory and was rebuilt by hand every round, so the "examples N/N"
in each PR had no authority in the tree.

Local ground-truth run: `matlab -batch "cd('dev/tests'); ci_gate"` — on the
user's PC preferably through the MATLAB MCP bridge (official
matlab-mcp-server; ~4-min call window: detach long runs as
`matlab -batch` jobs writing a result file, poll with `isfile`).

## Hard rules

- The gate is `ci_gate`, NEVER bare `runtests` (coverage identity §0a-2).
- No aliases; one name per thing; single point of access.
- One WP = one branch (`claude/wp-NN`) = one PR. Open a DRAFT PR right
  after the FIRST push (user decision 2026-08-19): `ci.yml` triggers on
  `pull_request` + push-to-main only, so the draft is what gives every
  dev push its ONE CI verdict (merge-preview context); mark it ready
  after both gates are green.
- The PR is opened ONLY after (1) CI green on the branch AND (2) the
  local MATLAB postflight is green; the postflight result is cited in
  the PR description. `git push --no-verify` is forbidden.
- The WP's OWN PR carries its records: the §1 ledger row, already marked
  DONE, and the dated round note in `HANDOVER.md`. There is NO separate
  closure PR (user decision 2026-11-08): it left `main` one records-only
  merge past every release tag, which toolboxInstaller reports as "ahead
  of the latest release". The row reaches `main` only through the merge
  that makes it true, so it is never early. DONE = merged to `main`,
  checks green; the row names the date the row is written (the session's
  day in the WP's branch, not the merge day, which is not known yet and
  is on the PR), the PR number and, when the PR bumps the version, the
  tag `vX.Y.Z` (release.yml puts it on the merge commit). It names no
  merge hash, which does not exist before the merge:
  `gh pr view NNN --json mergeCommit -q .mergeCommit.oid` finds it, exact
  for every merge style (not `git log --merges --grep "#NNN"`: the grep
  matches by prefix, `--merges` misses squash merges). A record that needs
  correcting after a merge is corrected in the next WP's PR, not alone.
- `main`'s ruleset requires STRICT status checks: green is not enough, the
  branch must also be up to date. After anything else lands, merge
  `origin/main` into the branch and let CI run again — "CI green but the
  merge button is dead" means BEHIND, not broken (2026-09-07).
- Every touched help block gets the Claude stamp (date, model) after
  `See also`. MATLAB coding standards apply
  (github.com/matlab/rules/blob/main/matlab-coding-standards.md).
- Measurement data and secrets NEVER enter the repo (public!).

## Token handling (Auth note, HANDOVER §5 WP-30 — binding)

The container authenticates with a fine-grained PAT (Contents RW,
Pull requests RW, Actions RW) stored as a chat-project file. Read it
from its file AT CALL TIME and pass it in a header; NEVER echo it into
command output, NEVER commit it, NEVER bake it into a remote URL that a
log could print. The user's own pushes use Git Credential Manager;
cloud sessions use the OAuth grant and need no token.

## References

Ellmer (2018) and Kvas (2019) are cited in `dev/references/README.md`
(canonical open-access records; `dev/references/` since L8). Full texts are
NOT redistributed in this public repository; for local sessions place
copies in `dev/references/local/` (gitignored). Equation numbers refer to Ellmer
unless said otherwise. `dev/references/BEST_PRACTICE_v4-1.md` Appendices
C/D are the operating runbook for container sessions.

## Repository layout since L8 (2026-09-02)

Everything development-related lives under `dev/`: `dev/HANDOVER.md` (SSOT),
`dev/tests/`, `dev/docs_src/` (audit + doc builders), `dev/references/`,
`dev/validation/`, `dev/SETUP_GUIDE.md`, `dev/HISTORY_v1.md`,
`dev/DEVELOPMENT.md` (the development chapters that used to sit in README).
One dev-related file STAYS at the root: this CLAUDE.md (the session entry
point by convention). `buildfile.m` was REMOVED in L9 (user decision
2026-09-02: an mltbx is not foreseen; the gate is `cd('dev/tests'); ci_gate`
and the audit is `python3 dev/docs_src/audit_tree.py` -- one way each). Tests manage their path
THEMSELVES, temporarily (run_all_lsqvce / ci_gate capture and restore it);
`lsqvce_setup.m` adds only root + examples.
