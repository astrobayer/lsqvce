# Mirror of lsqvce.util.log fired-inventory state (persistent FIRED) and the ORDER of calls in
# fit>runBatch (kinds B/C): warnOnce('fit:priorNotActive') at fit.m:274, then
# runBatchOnBlockPath -> inventory('reset') at fit.m:541 (or fit.m:289 on the x0 path),
# then the driver; R.meta.warnings = inventory() read by runMeta at the end.
FIRED = []
raised = []                      # what MATLAB's warning() actually shows
def warn(i): 
    if i not in FIRED: FIRED.append(i)
    raised.append(i)
def warnOnce(i):
    if i not in FIRED: warn(i)
def reset(): FIRED.clear()
def fit_kindB_with_late_prior():
    raised.clear()
    warnOnce('lsqvce:fit:priorNotActive')     # fit.m:274
    reset()                                   # fit.m:541
    meta = list(FIRED)                        # runMeta at the end (no other warning)
    return list(raised), meta
for run in (1, 2):
    shown, meta = fit_kindB_with_late_prior()
    print('run %d: warning raised=%s  R.meta.warnings=%s' % (run, shown, meta))
def fit_kindA_with_late_prior():
    raised.clear(); reset()                   # fit.m:172 resets FIRST
    warnOnce('lsqvce:fit:priorNotActive')     # fit.m:177
    return list(raised), list(FIRED)
print('kind A for contrast:', fit_kindA_with_late_prior())
# a kind-A run that fired the id leaves it in FIRED (nothing resets at the END of a run);
# the next kind-B run's warnOnce (before ITS reset) then sees it and stays silent.
raised.clear(); FIRED.clear()
fit_kindA_with_late_prior()
shown, meta = fit_kindB_with_late_prior()
print('kind A then kind B: kind-B warning raised=%s  kind-B R.meta.warnings=%s' % (shown, meta))
