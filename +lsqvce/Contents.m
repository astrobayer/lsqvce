% LSQVCE  Batch least squares: accumulation, variance components, and judgment.
% Version 6.3.0  10-Oct-2026
%
% lsqvce -- least-squares adjustment with variance component estimation.
%
% THE USER TIER (five names since R5-D8; HANDOVER §11.4). Start with the
% Manual (doc lsqvce): chapter 1 overview, chapter 3 the data you declare.
%
%   fit            - the one user verb: fit(data, opts) -> Results
%   layout         - the column layout fit will use (declare priors in it) (R5-D8)
%   Options        - the one settings object for fit (typed knobs with engine defaults)
%   Results        - rich result object: solution, diagnostics, provenance, report()
%   report         - verdict printer (levels 0-4); Results.report forwards here
%
% The engines and kernels below fit live in two sub-packages (R4-D1,
% 2026-09-01); they are reachable for experts and never needed by a user:
%   engine         - ENGINE tier: drivers, accumulators, loaders, VCE
%                    machinery, diagnostics    -> help lsqvce.engine
%   kernel         - KERNEL tier: pure numerical kernels (whitening,
%                    spectra, transforms)      -> help lsqvce.kernel
%   util           - shared utilities            -> help lsqvce.util
%
% Regenerated 2026-09-01 by docs_src/refactor_tiers.py (WP-64a).
