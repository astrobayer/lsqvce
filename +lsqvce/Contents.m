% LSQVCE  Batch least squares: accumulation, variance components, and judgment.
% Version 7.0.0  03-Nov-2026
%
% lsqvce -- least-squares adjustment with variance component estimation.
%
% THE USER TIER (v7.0.0: the declaration is a TaskDefinition; a data struct is the
% engine tier's, lsqvce.engine.fitDeclaration). New here? Run demo_tutorial and read
% the Manual's tutorial beside it. Start with the Manual (doc lsqvce) for how to use
% it, the API reference for the task you declare.
%
%   TaskDefinition - the declaration: flat add-calls compiled for fit (v7, WP-160)
%   SystemDefinition - one system of a task: locals, streams (Build), priors, Prepare
%   Noise          - the temporal law and level of a stream (white, arma, psd, ...)
%   DayContext     - what Prepare and Build know about their day
%   DayData        - the arcs one Build returns (addArc)
%   fit            - the one user verb: fit(task, opts) -> Results
%   layout         - the column layout fit will use (declare priors in it) (R5-D8)
%   Options        - the one settings object for fit (typed knobs with engine defaults)
%   Results        - rich result object: solution, diagnostics, provenance, report()
%   report         - verdict printer (levels 0-4); Results.report forwards here
%   solveContributions - a weighted combination of a fit's days (bootstrap) (WP-158)
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
