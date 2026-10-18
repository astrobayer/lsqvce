function dg = localDiagnostics(Nll, opt)
%LSQVCE.ENGINE.LOCALDIAGNOSTICS  Stability test of a local block before pre-elimination.
%
%   dg = LSQVCE.ENGINE.LOCALDIAGNOSTICS(Nll, opt) inspects the local normal matrix and
%   classifies its conditioning. Called by LSQVCE.ENGINE.BLOCKACCUMULATOR immediately before
%   the Schur complement, because a singular Nll silently corrupts the COMMON
%   block: pre-elimination forms
%
%       Ncc_red = Ncc - Nlc' * inv(Nll) * Nlc
%
%   so any junk in inv(Nll) propagates straight into the parameters you care
%   about. Diagnose before eliminating, never after.
%
%   INPUT
%     Nll  nl x nl  local normal matrix (symmetric)
%     opt  struct   .TruncTol   rank tolerance, relative to lambda_max
%                   .CondWarn   condition number -> 'ill-conditioned'
%                   .CondAbort  condition number -> 'rank-deficient'
%
%   OUTPUT  (WP-151: the spectrum is that of the JACOBI-EQUILIBRATED block
%           D Nll D, D = diag(Nll)^-1/2 -- the conditioning of the parameters
%           themselves, not of their units; the solve decides its rank on the
%           same matrix, LSQVCE.ENGINE.BLOCKACCUMULATOR)
%     dg.eig         nl x 1  eigenvalues of D Nll D, descending
%     dg.V           nl x nl its eigenvectors, matching the sorted eigenvalues
%     dg.condNll     spectral condition number lambda_max / lambda_min of D Nll D
%     dg.rankDef     number of eigenvalues below TruncTol*lambda_max
%     dg.spd         true if all eigenvalues are positive
%     dg.diagDom     min|diag| / max|diag|, a cheap scaling check
%     dg.nullDir     the weakest direction D v_min, unit norm -- the direction
%                    in local-parameter space that the data does not constrain
%     dg.status      'ok' | 'ill-conditioned' | 'rank-deficient'
%     dg.regUsed     'none' (filled in later by LSQVCE.ENGINE.APPLYREGULARISATION)
%     dg.lambdaUsed  0      (idem)
%
%   READING THE RESULT
%     dg.nullDir is the most actionable output. If it is concentrated on two
%     local parameters, those two are not separately estimable from this
%     system's data -- e.g. an accelerometer bias and a once-per-rev empirical
%     acceleration over a short arc. The fix is to remove one, not to
%     regularise: regularisation biases the common block through the Schur
%     complement.
%
%     dg.diagDom << 1 means the local parameters are badly scaled relative to
%     each other (e.g. bias in m/s^2 alongside a dimensionless scale factor).
%     A large condNll caused purely by scaling is an artefact -- condNll reads
%     the equilibrated block, so it no longer shows one (WP-151) -- but
%     Tikhonov with L = I penalises the two parameters unequally: rescale the
%     columns of B, or give RegMatrix the scales.
%
%   TYPICAL THRESHOLDS
%     CondWarn  = 1e8   double precision starts to bite past this
%     CondAbort = 1e14  effectively singular; 'auto' regularisation goes TSVD
%
%   EXAMPLE
%     opt = struct('TruncTol',1e-10, 'CondWarn',1e8, 'CondAbort',1e14);
%     dg  = lsqvce.engine.localDiagnostics(Nll, opt);
%     if ~strcmp(dg.status,'ok')
%         lsqvce.util.log.msg('cond = %.2e, weakest direction:\n', dg.condNll);
%         disp(dg.nullDir.')
%     end
%
%   ACCURACY
%     No reference value is pinned: the spectral outputs come from one
%     symmetric eig of D (Nll + Nll')/2 D and serve as a classification. A
%     local column in seconds squared (B = [1, t^2], t to 2e4 s, raw cond
%     ~1e17) reads 'ok' (test_lsqAccum (testRankDecisionIsUnitFree)). Pinned
%     by behaviour: near-collinear local columns (collinearEps 1e-9) read
%     condNll > 1e8 and a status other than 'ok'; at 1e-12 nullDir sits on
%     the two collinear columns with opposite signs; a zero 3 x 3 block
%     reads rankDef 3 (test_reliability (testDiagnosticsFlagIllConditioning,
%     testNullDirectionIdentifiesTheCulprit,
%     testAZeroLocalBlockIsWhollyRankDeficient)). WP-139.
%
%   ERRORS
%     lsqvce:localDiagnostics:notFinite  a local block carries NaN or Inf
%
%   See also LSQVCE.ENGINE.APPLYREGULARISATION, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-151 (review 12.84, S1): the spectrum, condNll, rankDef,
%   status and nullDir read the Jacobi-equilibrated block, so a block in mixed
%   units is not classified ill-conditioned (and 'auto' does not regularise it).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3, lane 1): rankDef
%   counts eigenvalues at or below TruncTol * max(lambda_max, 0) -- an all-zero
%   block read 0 deficient directions while its status said rank-deficient. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(2,2);
    validateattributes(Nll, {'double'}, {'2d','square','nonempty'}, mfilename, 'Nll', 1);

    if ~all(isfinite(Nll(:)))
        error('lsqvce:localDiagnostics:notFinite', ...
            ['The local normal matrix contains NaN or Inf. The design matrix ' ...
             'or the whitening produced non-finite values -- check for a zero ' ...
             'sigma, an empty arc, or a gap in the observations.']);
    end

    tolT = lsqvce.util.getfielddef(opt, 'TruncTol',  1e-10);
    cW   = lsqvce.util.getfielddef(opt, 'CondWarn',  1e8);
    cA   = lsqvce.util.getfielddef(opt, 'CondAbort', 1e14);

    Nll = (Nll + Nll.')/2;
    dd  = abs(diag(Nll));
    % WP-151 (S1): the spectrum of the Jacobi-equilibrated block D Nll D -- a
    % condition number "caused purely by scaling is an artefact" (above), and
    % 'auto' regularised a well-determined t^2 local (t in seconds) on it
    d   = lsqvce.engine.BlockAccumulator.jacobiScale(Nll);
    Nll = Nll .* (d * d.');
    [V, ev] = eig(Nll, 'vector');              % symmetric -> V orthonormal
    ev      = real(ev);
    [ev, k] = sort(ev, 'descend');
    V       = real(V(:, k));

    dg          = struct();
    dg.eig      = ev;
    dg.V        = V;
    dg.condNll  = ev(1) / max(abs(ev(end)), realmin);
    dg.rankDef  = sum(ev <= tolT*max(ev(1), 0));   % WP-135: an all-zero block is nl deficient (it read 0)
    dg.spd      = all(ev > 0);
    dg.nullDir  = d .* V(:, end);              % back in the parameters' units
    dg.nullDir  = dg.nullDir / max(norm(dg.nullDir), realmin);

    dg.diagDom  = min(dd) / max(max(dd), realmin);

    dg.status = 'ok';
    if dg.condNll > cW,              dg.status = 'ill-conditioned'; end
    if dg.condNll > cA || ~dg.spd,   dg.status = 'rank-deficient';  end

    dg.regUsed    = 'none';
    dg.lambdaUsed = 0;
end
