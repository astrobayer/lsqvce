% lsqvce.engine -- the ENGINE tier of lsqvce (R4-D1, WP-64a).
%
% Users reach these through lsqvce.fit / lsqvce.Options / lsqvce.Results;
% experts call them directly. Roles from docs_src/tiers.csv, the ONE tier
% authority (audit check 18).
%
%   adjustNonlinear    - Gauss-Newton / Levenberg as a policy on the block path (R5 M2)
%   batchSystem        - a batch source as one local-free streaming system (R5 M2)
%   adjustRobust       - IRLS (Huber/IGG3/bisquare) as a policy on the block path (R5 M2)
%   adjustBlock          - one GLS solve over the streaming block-arrow pipeline
%   adjustVCEjoint       - Ellmer ten-step joint co-estimation with keyed arc factors
%   applyRegularisation  - regularisation of the reduced system
%   arcVarianceFactors   - per-arc variance factors
%   BlockAccumulator     - block-arrow accumulation with per-system Schur pre-elimination
%   BlockCofactor        - the cofactor of the full parameter set, kept as blocks (WP-86)
%   checkPSD             - PSD sanity check
%   collectResiduals     - residual collection per stream
%   CommonPrior          - prior on the common block (regularisation as observations)
%   contribution         - information contribution per group
%   covComponent         - covariance component estimation
%   covModel             - covariance spec -> whitening operators (per arc)
%   datum                - datum defect detection and handling
%   declaration        - validate a kind-A data struct and build its loaders (R5-D8)
%   gcvLambda            - GCV choice of the regularisation weight
%   jointModelOpt        - joint-mode model -> loader UserOpt translation
%   localDiagnostics     - diagnostics of the local blocks
%   localColumns       - which (row, day, arc) owns each column of B (R5-D8)
%   makeLoader           - day-builder -> Loader (built INSIDE fit from the data struct, R4-D4)
%   mapreduce            - serial / parfor / spmd task runner with ordered reduce
%   NormalEqSystem       - normal-equation backend
%   outliers             - outlier listing from diagnostics
%   psdVCE               - PSD amplitude + arc-factor VCE (Ellmer workhorse)
%   QRSystem             - QR backend
%   ReducedContribution  - reduced-normal contribution bookkeeping
%   redundancy           - redundancy numbers (data pass)
%   residuals            - second data pass: residuals and quadratic forms
%   resolveNc            - common-parameter count derived from the systems
%   synthLoader        - synthetic block-arrow loader (engine-tier fixture; leaves in M4)
%   sharedProbe          - shared Hutchinson probe generator (RNG contract)
%   streamContribution   - contribution per stream
%   traceEstimate        - stochastic trace estimate
%   vce                  - Helmert group variance components (exact | hutchinson)
%   verdict              - the verdict: diagnostics, snooping, systematics, F test, ellipse (WP-86)
%   runMeta              - immutable provenance record for Results.meta
