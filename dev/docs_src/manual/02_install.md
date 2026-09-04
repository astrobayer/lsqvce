# Installation

## Requirements

MATLAB R2024a or later (developed and gated on R2026a). The Parallel Computing Toolbox
is optional (`Parallel.Mode="parfor"|"spmd"`, `Parallel.DistributeDays`); without it the
toolbox runs serially and says so once per run.

## Install

```matlab
cd <toolbox root>          % the folder holding +lsqvce
lsqvce_setup               % or lsqvce_setup('install')
```

This adds exactly two path entries (the root with the `lsqvce.*` packages, and
`examples/`), builds the Help browser index so that `doc lsqvce` opens this manual, saves
the path, and prints a check report. Not added by design: `dev/tests` (the suites manage
their own path), `template/` (you copy it into your project), `dev/`.

## Check

```matlab
lsqvce_setup('check')      % version, path entries, doc index, whether lsqvce.fit resolves
```

## Uninstall

```matlab
lsqvce_setup('uninstall')  % removes the two path entries, saves the path, deletes the
                           % Help browser database it built; the folder itself is untouched
```

## Layout of the toolbox folder

| folder / file | what |
|---|---|
| `+lsqvce/` | the package: user tier (`fit`, `Options`, `Results`, `report`, `layout`), `+engine`, `+kernel`, `+util` |
| `examples/` | the demos (`demo_workflow`, `demo_fit`, `demo_grace`, `demo_gravity`, `demo_gravfire`, ...) |
| `template/` | `gravfire_data.m` (a GRACE bundle to the data struct), `lsqvce_template.m` |
| `html/` | this manual (`lsqvce_manual.html`) and one page per public unit (`doc lsqvce.<name>`) |
| `lsqvce_theory.pdf`, `lsqvce_manual.pdf` | the theory and this manual as PDF |
| `RELEASE_NOTES.md` | what changed, with the migration table |
| `dev/` | the development record (`HANDOVER.md`), the test suites, the documentation sources |
