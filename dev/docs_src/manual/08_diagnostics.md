# Diagnostics and the verdict

Every solve carries per-observation diagnostics in `R.diag` (residuals in the whitened
frame `sqrt(p).*v`, redundancy numbers, standardised residuals, the Baarda / Pope test
with `Diagnostics.Alpha` and Bonferroni, minimal detectable biases, reliability, the
a-priori backstop). `R.report(Level=0..4)` reduces them to one verdict and its evidence;
`R.summary()` and `R.convergence()` read the pass record of an iterative run;
`R.residuals()` re-reads the arcs for the whitened residuals per stream; `R.contribution()`
gives the information share of every stream in the common parameters;
`lsqvce.engine.snoop` iterates data snooping with a clean re-solve on the block path;
`R.compare(Rred)` is the nested-model F test between two fits. For correlated batches the
per-observation quantities live in the whitened frame; only their sums are
frame-invariant.
