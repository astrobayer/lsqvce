# Common rules for every reviewer lane (lsqvce review R7, 2026-10-09)

Repository (read-only for you): /home/claude/lsqvce  -- main@b25bb11, toolbox version 6.1.3 (+lsqvce/Contents.m).
Skill dir: /root/.claude/skills/synced/4108e6bc-08ac-477c-9f82-10f57d4b005c_6628b5d9-0ff0-4a82-8155-817a1ef1a1c7/scientific-toolbox-discipline
  (references/review_roles.md = output format + common rules; references/coding.md, validation.md, documentation.md).

TIER: B. There is NO MATLAB here. Python 3 with numpy 2.5 and scipy 1.18 IS available. A finding is REPRODUCED
only if you computed it (a Python mirror of the arithmetic, a grep whose output proves it, running the repo's Python
audit). Everything derived from reading is INFERRED. Never mix the two. Python mirrors go in your own scratch dir
(/tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/mirrors/<lane>/), never in the repo.

DO NOT modify any file under /home/claude/lsqvce. Do not git commit/push. Findings only.

Single point of truth: dev/HANDOVER.md (7643 lines). Read at least: §0 (lines 371-434), §12.2 (R5 confirms as RIGHT),
§12.6 decisions, §12.14 (R6 review, findings E-01..E-13), and the last rounds §12.29-§12.34 (lines ~7072-7643).
A recorded decision is NOT a finding; CONTRADICTING one is (say so). Do not re-propose closed findings unless the
code shows they regressed. Also: dev/DEVELOPMENT.md, dev/docs_src/error_ids.csv, dev/tests/category_exemptions.csv.

Already found by the coordinator (do not re-report; you may extend):
 C1 adjustVCEjoint keeps the claim WP-110 retired from adjustRobust: help lines 175-177 "(beyond the breakdown point")
    split over two lines so audit check 28's single-line regex misses it, and the runtime message at line ~649-653
    ("more than half the data is discrepant, which is beyond the breakdown point") -- check 28 does not scan error
    message bodies. Its guard also compares live EPOCHS with nc+1 PARAMETERS.
 C2 lsqvce:adjustVCEjoint:badSigma2_0 and :arcMismatch are raised and catalogued in error_ids.csv but absent from the
    function's ERRORS block and fired by no test.

REPORT: write ONE markdown file to
 /tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/reviews/<lane>.md
in exactly the shape of review_roles.md "Output format": Findings most severe first, each with
 title — REPRODUCED|INFERRED; What is wrong; Why; Observable symptom (falsifiable); Evidence (file:line, command+output);
 Severity: S1 (silently wrong numbers) / S2 (wrong under reachable conditions, or a refusal/guard that misfires) /
 S3 (contract/doc/test gap) / S4 (hygiene);
then "Checked and sound" (one line per area, how verified), then "Not reviewed".
Quality over quantity: at most ~15 findings, each verified against the CURRENT code (quote the line). No speculation
presented as fact. Your final message back: a 10-line summary of the top findings with severities and the report path.
