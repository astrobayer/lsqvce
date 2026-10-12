# Installation

## Requirements

MATLAB R2024a or later. The Parallel Computing Toolbox is optional (`Parallel.Mode =
"parfor" | "spmd"`); without it the toolbox runs serially and says so once per run.

## Install

```matlab
cd <toolbox root>          % the folder holding +lsqvce
lsqvce_setup               % or lsqvce_setup('install')
```

This adds exactly two path entries — the root with the `lsqvce.*` packages, and
`examples/` — builds the Help browser index so that `doc lsqvce` opens this Manual, saves
the path, and prints a check report. `template/` is not added: you copy a template into
your own project.

## Check

```matlab
lsqvce_setup('check')      % version, path entries, doc index, whether lsqvce.fit resolves
```

## Uninstall

```matlab
lsqvce_setup('uninstall')  % removes the two path entries and the Help browser index;
                           % the folder itself is untouched
```

## What is in the folder

| folder / file | content |
|---|---|
| `+lsqvce/` | the package: the user tier (`fit`, `Options`, `Results`, `report`, `layout`) and the `engine`, `kernel` and `util` tiers below it |
| `examples/` | runnable examples (chapter *The examples*) |
| `template/` | starting points for your own case: `kbr_variational.m` (GRACE-FO Level-1B range-rates and orbits), `gravfire_data.m` (a GRAVFIRE bundle to the `data` struct), `lsqvce_template.m` |
| `html/` | this Manual, the API reference and one page per public unit |
| `lsqvce_manual.pdf`, `lsqvce_api.pdf`, `lsqvce_theory.pdf` | the three documents as PDF |
| `RELEASE_NOTES.md` | what changed in each version |
| `dev/` | development: tests, documentation sources, the development record |
