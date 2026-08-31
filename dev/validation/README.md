# Python validation mirror

These scripts are the numerical ground-truth harness used to validate the
MATLAB toolbox before transcription (MATLAB is not available in the build
environment). Each brute-forces the quantity the corresponding MATLAB code
computes, and asserts agreement.

    python3 verify.py      # v2 core: 23 checks
    python3 verify2.py     # v2 precision + reliability: 24 checks
    python3 v3verify.py    # v3 features standalone: 18 checks
    python3 v3integ.py     # v3 integrated through the accumulator: 9 checks
    python3 v4verify.py    # v4 parallel decomposition + iterateNoiseModel: 8 checks
    python3 v4boundary.py  # dense-path boundary-bias MEASUREMENT (2.2)
    python3 v5verify.py    # v4.1 day-distributed accumulation: 6 checks
    python3 v6verify.py    # v4.2 VAR path + spmd tree combine: 8 checks

v3core.py holds the reference implementations imported by v3verify / v3integ.

v4verify.py proves the decompositions the MATLAB parallel paths rely on: arc-
and system-level partition sums reconstruct the serial normals exactly;
task-partitioned redundancy equals the fused serial loop in any task order;
iterateNoiseModel converges from a white start AND recovers the true model ACV; the
residual-ACV leverage bias and the quick-look bias `1 - sqrt(r/m)` are
measured. v4boundary.py is not a pass/fail test but the measurement behind the
"BatchSize >= 4x correlation length" rule: no cliff, no bias, `b^(-1/2)` decay
of the boundary perturbation, efficiency loss <1% at b = 4L.

v5verify.py proves the day-distributed accumulation (v4.1), in particular the
QR path: compressed (dim+1)-row triangular factors are lossless in N, n and
l'l; per-day factors merged by stack + re-QR equal the fused serial
accumulation to ~5e-13 in every derived quantity, in any merge order; the
residual-row mass (what keeps sigma0 honest) survives the merge -- and the
test verifies the guard has teeth by showing what happens when the row is
dropped; the full day-distributed pipeline matches a monolithic dense
reference for both backends.

Working rule for all future development: validate a new feature here against a
brute-force reference BEFORE writing the MATLAB, and re-run the full set after.
This has caught every real bug in the project.

## Note on the two statistical tests

`verify2.py`'s "formal == empirical sigma" and `v3integ.py`'s "VCE recovers
variance components" are STATISTICAL checks over random noise realisations.
Their tolerances are calibrated to the sampling scatter, so on an occasional RNG
draw one may land just over threshold (e.g. 23/24 instead of 24/24). Re-run and
it passes. This is expected finite-sample behaviour, not a code fault — the mean
behaviour (the actual claim: formal errors calibrated, components recovered) is
correct every run. All other checks are deterministic and pass every time.

v6verify.py proves the VAR path before transcription: the companion-form
cross-ACV satisfies Yule-Walker exactly and matches simulation; the Whittle
recursion recovers the generating (A, Sigma) from the exact ACV and overfit
orders vanish; the matrix-FIR whitener with exact block startup equals the
dense block-Cholesky inverse ROW FOR ROW (~1e-10) and satisfies the GLS
identities; batched whitening with p-epoch state equals whole-arc; and the QR
mergeRaw combine is grouping-independent, which is what makes the spmd tree
reduce of day partials safe.
