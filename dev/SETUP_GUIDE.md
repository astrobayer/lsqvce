# lsqvce — Setup Guide: public GitHub + hosted CI + the cloud-session loop (D8)

*For a first-time GitHub/CI user. Follow it top to bottom.
Design and status live in `HANDOVER.md` (WP-30, decision D8, corrected);
this file is HOW-TO only. Part 6 shows how to reuse the pattern for any
other project.*

*Created 2026-08-17 (CEST, UTC+02:00) with assistance from Claude
(Anthropic), model Claude Fable 5. Replaces the earlier same-day guides
(Route C, D6, first D8 draft).*

---

## The picture (read once, 2 minutes)

**The loop:**

    THINK      in Claude chat / Projects (or Cowork): discuss, design,
               pin the work-package prompt into HANDOVER.md.
    PROMPT     hand that prompt to Claude in chat or Cowork — the
               session's own container clones the repo, opens branch
               claude/wp-NN, develops, and pushes. No pull request
               yet. (Token-free alternative: a cloud session.)
    CI         every push gets the CI verdict (checks: audit + matlab,
               real MATLAB on GitHub's Linux machines). The session
               polls the checks itself, reads red logs, fixes. To green.
    POSTFLIGHT CI green? At the PC, Claude checks out the branch and
               runs the SAME gate in YOUR MATLAB through the MCP
               bridge — second platform: Windows, your R2026a.
    PR         both green? Only now the pull request is opened, with
               the postflight result in its description. You read the
               diffs and logs, comment, then merge. Merging is yours.
    RECORD     the same PR updates HANDOVER.md's ledger:
               DONE = merged to main, checks green, date + PR + hash.

Five actors: **GitHub** (code, history, CI, PRs — public), **Claude
cloud sessions** (the developer), **CI MATLAB** (the judge, free
because the repo is public), **Claude chat/Projects/Cowork** (where you
think and prompt), **you** (the only one who merges). Your PC is
optional; your local MATLAB stays your personal workbench — and, via
the MCP bridge, an optional fast test path.

**Public means public.** The free MATLAB on GitHub's runners exists
only for public repositories — that is the deal. Consequences: anyone
can read the code; measurement data and secrets NEVER enter the repo
(the `.gitignore` enforces the data part); and you should decide on a
LICENSE file (until you do, "public, all rights reserved" applies —
people can look, not legally reuse).

**What you need:** a GitHub account, your Claude account (Pro or
higher), MATLAB R2026a on the PC (for yourself and the bridge).

---

## Part 1 — Base tools on your PC (~20 min)

**1.1 Git:** install from https://git-scm.com/download/win (accept all
defaults). Check in a **new** PowerShell: `git --version`.

**1.2 MATLAB on the system PATH:** Windows key → `environment
variables` → *Edit the system environment variables* → *Environment
Variables…* → System variables → `Path` → *Edit* → *New* →

    C:\Program Files\MATLAB\R2026a\bin

Check in a new PowerShell: `matlab -batch "disp(version)"` prints an
R2026a version line (first start takes a minute).

**1.3 Claude Desktop app:** update if needed (claude.ai/download), sign
in. The *Chat* tab is where you think and prompt; the *Code* tab
(Claude Code with a GUI) is optional for local sessions.

**1.4 The MATLAB MCP bridge — the postflight stage:** this connects
Claude Desktop to your local MATLAB. In the loop it is the mandatory
POSTFLIGHT check: after CI is green, the branch is verified once in
your own MATLAB before a pull request may be opened. (It also works as
an optional pre-flight during development.)

- In Claude Desktop: *Settings* → *Extensions* → *Browse extensions* →
  install **Filesystem** (by Anthropic) and allow it the folder
  `C:\dev\lsqvce`.
- Download `matlab-mcp-server.mcpb` from the Latest Release at
  https://github.com/matlab/matlab-mcp-server (official MathWorks,
  free) and double-click it → *Install* in Claude Desktop.
- Restart the Desktop app; the tools icon should list MATLAB tools
  (run code, check code, start/quit MATLAB).

> **Checkpoint 1:** `git --version` and `matlab -batch "disp(version)"`
> work; (optional) in the Chat tab, "start MATLAB and run disp(2+2)"
> succeeds through the bridge.

---

## Part 2 — The repository (~10 min)

**2.1** github.com → **+** → *New repository* → name `lsqvce` →
visibility **Public** → no "initialize" boxes → create.

**2.2** Unzip the current `lsqvce.zip` to `C:\dev\lsqvce` so that this
folder directly contains `+lsqvce\`, `tests\`, `HANDOVER.md`, and this
guide.

**2.3** First upload (PowerShell; git opens a GitHub sign-in once):

    git config --global user.name  "Matthias ..."
    git config --global user.email "you@example.org"

    cd C:\dev\lsqvce
    git init
    git add .
    git commit -m "R2 baseline"
    git tag r2-baseline
    git branch -M main
    git remote add origin https://github.com/<YOU>/lsqvce.git
    git push -u origin main --tags

> **Checkpoint 2:** the GitHub page shows your files and the tag
> `r2-baseline`.

---

## Part 3 — The sessions' GitHub credential (~10 min)

The developer in this loop is the chat/Cowork session's own container —
your proven Appendix-D pipeline (`references/BEST_PRACTICE_v4-1.md`).
A container has no GitHub identity, so it needs the token you already
know:

- GitHub → *Settings* → *Developer settings* → *Fine-grained tokens* →
  *Generate new token*: Repository access = **Only** `lsqvce`;
  Permissions: **Contents: Read and write**, **Pull requests: Read and
  write**, **Actions: Read and write** (Actions is what lets a session
  cancel and re-run a hung CI job — a missing scope is discovered
  mid-incident, which is too late). Choose an expiry you can live with.
- Store it as a **project file** in the Claude Project (e.g.
  `github_token.txt`); sessions read it at call time.
- Handling rules (CLAUDE.md restates them): never echoed into command
  output, never committed, never baked into a remote URL that a log
  could print.

Your OWN pushes still just sign in once via the browser (Git Credential
Manager); public clones need no authentication. *Token-free
alternative:* Claude Code cloud sessions via the one-time GitHub grant
at claude.ai/code — handy from the phone, less scriptable.

> **Checkpoint 3:** in a fresh chat session, one probe API call with
> the token returns your repo's metadata.

---

## Part 4 — The pipeline files (one Claude round)

**4.1** Start a **new** Claude chat, paste the §9.0 preamble plus the
WP-30 section from `HANDOVER.md`, attach the current zip. That round
delivers `.github/workflows/ci.yml` (checks `audit` + `matlab`, MATLAB
R2026a on GitHub's machines), `tests/ci_gate.m`, `CLAUDE.md`,
`.gitignore`, `.gitattributes`, `references/README.md` (citations; the
thesis texts stay OUT of the public repo — drop local copies into the
gitignored `references/local/` for repo-based sessions) — and asks you
the LICENSE question.

**4.2** Commit and push; the **Actions** tab must show `audit` and
`matlab` green on `main`. The very first `matlab` run installs MATLAB
on the runner and takes several minutes — that is normal.

**4.3 The red test:** push the deliberately failing commit the WP-30
round prepares, on a branch: the checks must go RED and the PR must
show as blocked. Revert.

**4.4 Branch protection:** *Settings* → *Branches* → rule for `main` →
require the status checks `audit` and `matlab` before merging.

> **Checkpoint 4:** green on `main`; a red branch is visibly
> unmergeable.

---

## Part 5 — Daily use: the loop itself

1. **Think** in Claude chat or the Project — exactly like the sessions
   that produced `HANDOVER.md`. Outcome: a pinned work-package prompt.
2. **Prompt:** hand it to Claude in chat or Cowork. The session's
   container fetches, checks branch/PR state, opens `claude/wp-NN`,
   pre-validates in Python, implements, and pushes. No pull request
   yet. (From the phone without a container: a cloud session.)
3. **CI to green:** the CI judges every push; the session polls the
   check-runs itself at suite cadence, reads red logs, and fixes — you
   step in only where judgment is needed.
4. **Postflight:** at the PC, in chat with the bridge active:
   "check out branch wp-33 and run ci_gate through the MATLAB bridge."
   Same gate, your Windows R2026a. If it disagrees with the CI, that
   disagreement is a finding — back to step 3.
5. **PR:** both green → Claude opens the pull request, citing the
   postflight result in the description. You read *Files changed* and
   the log artifact, request changes as comments, then **Merge**.
6. **Record:** the PR contains the HANDOVER ledger update (CLAUDE.md
   obliges it): DONE + date + PR number + short hash. Missing? Ask for
   it before merging.

Steps 1–3 work from the phone; postflight waits until you are at the
PC — a CI-green branch simply queues there.

---

## Part 6 — Transfer the pattern to any other project

The pattern is four files plus two clicks, independent of MATLAB:

1. A single **gate entry point** the CI calls (here `tests/ci_gate.m`;
   elsewhere maybe `make gate` or `pytest -x`): one command, green or
   fail, checks its own coverage.
2. **`ci.yml`** — copy from this repo; swap the gate command and the
   toolchain setup step (setup-matlab → setup-python, etc.).
3. **`CLAUDE.md`** — pointer to THAT project's single source of truth,
   its commands, its hard rules.
4. A **`.gitignore`** that keeps data and generated files out.
Then: make the repo visible to your Claude GitHub access, protect
`main` — and the same THINK → PROMPT → REVISE → REVIEW → RECORD loop
applies. Local toolchains plug in the same way MATLAB did: through an
MCP server, if one exists for them.

---

## If something goes wrong

- **The `matlab` check fails immediately with a licensing message** —
  is the repository still Public? The free MATLAB exists only for
  public repos.
- **Cloud session cannot see the repo** — redo Part 3 (GitHub access
  grant covers what your GitHub account can see).
- **The bridge does nothing** — Desktop app fully restarted after the
  install? Filesystem extension enabled for `C:\dev\lsqvce`? MATLAB on
  the system PATH (1.2)?
- **`git`/`matlab` "is not recognized"** — old terminal; open a new
  one, or redo 1.2.
- **A PR has no ledger update** — ask Claude in the PR/chat to add it
  before you merge.

## Four rules that keep this safe

1. **Public means public:** no measurement data, no credentials, no
   personal information in the repo — ever. When in doubt, it stays
   out.
2. Decide the LICENSE consciously; until then it is "look, don't
   reuse".
3. Any change under `.github/workflows/` is read by YOU before it is
   merged.
4. Merging is yours alone.
