function S = systematic(D, opt)
%LSQVCE.ENGINE.SYSTEMATIC  Detect unmodelled systematic errors in the residuals.
%
%   S = LSQVCE.ENGINE.SYSTEMATIC(D) runs three complementary tests on the normalised
%   residuals from lsqvce.engine.diagnose. This closes the single largest blind spot
%   in the rest of the toolbox.
%
%   WHY THIS EXISTS
%   ---------------
%   A smooth unmodelled signal -- a thermal drift, an unmodelled bias that
%   varies with elevation, a once-per-revolution term you forgot -- produces NO
%   single large residual. Data snooping therefore reports a perfectly clean
%   result while the parameters are quietly corrupted. Meanwhile the global
%   chi-square test rejects, because the residuals ARE too large in aggregate.
%
%   That combination -- inflated variance factor, zero outliers -- is the
%   signature of a systematic error, and until now the toolbox could only tell
%   you to go and plot something. Demonstrated on a sinusoid of amplitude 0.8
%   sigma: s0^2 = 1.33 (inflated), Bonferroni snooping flagged 0 observations,
%   and the three tests below returned p = 5e-5, p < 1e-16 and p = 3e-9.
%
%   THE THREE TESTS
%   ---------------
%   (a) RUNS TEST on sign(w). Distribution-free. Counts sign changes: a drift
%       produces too few. Weak on smooth signals, strong on steps.
%
%   (b) LJUNG-BOX on w, lags 1..h. Detects autocorrelation of any shape.
%       Powerful, but see the caveat below.
%
%   (c) BINNED MEANS of w against each supplied covariate. In bin b the mean
%       of w has standard deviation 1/sqrt(n_b) under H0, so z_b = mean_b *
%       sqrt(n_b). Bonferroni over the bins. This is effectively a score test
%       for an omitted parameter, and it is the WORKHORSE: in testing it
%       reached 84-100% power where the runs test managed 1-9%.
%
%   (a) and (b) require the residuals to be ORDERED along the systematic --
%   they see the sequence you hand them. (c) does not: it works on any
%   covariate directly. Pass covariates and let (c) do the work.
%
%   CAVEAT ON LJUNG-BOX -- and it is an honest one.
%   Least-squares residuals are NOT independent, even when the observation
%   errors are: cov(v) = sigma0^2 * Qvv, and Qvv is not diagonal. Normalising
%   to w removes the variance but not the correlation. Ljung-Box therefore
%   OVER-REJECTS, and the effect grows with u/n. Measured false-alarm rate at a
%   nominal 1%:
%       u/n = 0.060  ->  1.6%
%       u/n = 0.012  ->  1.3%
%       u/n = 0.003  ->  0.6%   (within noise of nominal)
%   Geodetic problems live at u/n << 0.01, where it is harmless. The toolbox
%   warns when u/n > 0.05, and it leads with the binned test, which is exact.
%
%   NAME-VALUE
%     'covariates'  struct: each field is an n x 1 vector to bin against.
%                   e.g. struct('time', t, 'elev', e, 'azim', az)
%                   THIS IS THE IMPORTANT ONE. Supply everything you have.
%     'nbins'       bins per covariate (default 10)
%     'lags'        Ljung-Box lags h (default 10)
%     'alpha'       significance (default 0.01)
%     'verbose'     print (default true)
%
%   OUTPUT S
%     .runs      struct(z, p, nRuns, expected)
%     .ljungBox  struct(Q, p, h, reliable)
%     .bins      struct array, one per covariate:
%                (name, z, p, worstBin, worstZ, edges, counts, means)
%     .detected  true if any test rejects at alpha
%     .worst     name of the most significant finding
%
%   INPUT
%     D         1x1 struct  diagnostics from LSQVCE.ENGINE.DIAGNOSE
%     varargin  Name/Value: 'covariate' (a vector to test against), 'alpha', ...
%
%   OUTPUT
%     S  1x1 struct  .runs, .durbinWatson, .groupOffsets, .verdict -- the
%                    structure tests that distinguish a MODEL defect from
%                    scattered blunders
%
%   EXAMPLE
%       % The failure that snooping CANNOT see. A smooth drift produces no
%       % single large residual: the variance factor inflates, the flag list
%       % stays clean, and THAT COMBINATION is the signature.
%       [~, D] = lsqvce.engine.adjust(u, src, 4, 'verbose', false);
%
%       t = (1:n).';
%       S = lsqvce.engine.systematic(D, 'covariates', ...
%                              struct('time', t, 'elev', elev));
%
%       % S.detected     did anything fire?
%       % S.ljungBox.p   autocorrelation in the residual sequence
%       % S.runs.p       distribution-free sign test
%       % S.bins         one entry per covariate: .name, .p, .worstBin
%
%       % The COVARIATE ranking lives in S.bins -- S.worst ranks every test,
%       % and a drift in time necessarily also trips Ljung-Box.
%       [~, iw] = min([S.bins.p]);
%       lsqvce.util.log.msg('worst covariate: %s\n', S.bins(iw).name);
%
%   ERRORS
%     lsqvce:systematic:noResiduals   the diagnostics carry no residuals
%     lsqvce:systematic:badCovariate  the covariate is not a vector matching the
%                                     residual length
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.FTEST.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (b) b2: typed `arguments` signature (canonical
% UpperCamelCase; legacy casings keep matching -- MATLAB NV
% semantics; enum VALUES case-sensitive lowercase).
arguments
    D
    opt.Covariates (1,1) struct = struct()
    opt.NBins (1,1) double {mustBeGreaterThanOrEqual(opt.NBins,2)} = 10
    opt.Lags (1,1) double {mustBeGreaterThanOrEqual(opt.Lags,1)} = 10
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 0.01
    opt.Verbose (1,1) logical = true
end

o = struct('covariates', opt.Covariates, 'nbins', opt.NBins, 'lags', opt.Lags, 'alpha', opt.Alpha, 'verbose', opt.Verbose);

if ~isfield(D, 'w')
    error('lsqvce:systematic:noResiduals', ...
        ['The diagnostics struct has no residual vector. Re-run ' ...
         'lsqvce.engine.diagnose with ''store'',''all'' -- the systematic tests ' ...
         'need every residual, not just the flagged ones.']);
end

w = D.w(:);
n = numel(w);
u = D.u;

% Only observations the network actually controls carry information about a
% systematic. An uncontrolled observation has w = 0 by construction, and
% including those zeros would dilute every test toward "no signal".
if isfield(D, 'rNum')
    ctrl = D.rNum(:) > D.rMin;
else
    ctrl = true(n, 1);
end
wc = w(ctrl);
nc = numel(wc);

S = struct();
S.n = n; S.nControlled = nc; S.alpha = o.alpha;

if nc < 20
    lsqvce.util.log.warn('lsqvce:systematic:tooFew', ...
        ['Only %d controlled observations. The systematic-error tests are ' ...
         'asymptotic and are meaningless at this size.'], nc);
    S.detected = false;  S.worst = 'n/a';
    S.runs = struct('z', NaN, 'p', NaN);
    S.ljungBox = struct('Q', NaN, 'p', NaN, 'reliable', false);
    S.bins = struct([]);
    return;
end

%% ------------------------------------------------------------- (a) runs test
sg = sign(wc);
sg = sg(sg ~= 0);
ns = numel(sg);
npos = sum(sg > 0);
nneg = ns - npos;

if npos == 0 || nneg == 0 || ns < 3
    S.runs = struct('z', NaN, 'p', NaN, 'nRuns', NaN, 'expected', NaN);
else
    nRuns = 1 + sum(sg(2:end) ~= sg(1:end-1));
    muR   = 2 * npos * nneg / ns + 1;
    varR  = 2*npos*nneg * (2*npos*nneg - ns) / (ns^2 * (ns - 1));
    if varR <= 0
        S.runs = struct('z', NaN, 'p', NaN, 'nRuns', nRuns, 'expected', muR);
    else
        z = (nRuns - muR) / sqrt(varR);
        p = 2 * (1 - lsqvce.engine.stat.normCdf(abs(z)));
        S.runs = struct('z', z, 'p', p, 'nRuns', nRuns, 'expected', muR);
    end
end

%% -------------------------------------------------------------- (b) Ljung-Box
h = min(o.lags, floor(nc / 5));
if h < 1
    S.ljungBox = struct('Q', NaN, 'p', NaN, 'h', 0, 'reliable', false);
else
    wz = wc - mean(wc);
    c0 = (wz' * wz) / nc;
    if c0 <= 0
        S.ljungBox = struct('Q', NaN, 'p', NaN, 'h', h, 'reliable', false);
    else
        Q = 0;
        for k = 1:h
            ck = (wz(1:end-k)' * wz(k+1:end)) / nc;
            rk = ck / c0;
            Q  = Q + rk^2 / (nc - k);
        end
        Q = nc * (nc + 2) * Q;
        p = 1 - lsqvce.engine.stat.chi2Cdf(Q, h);

        % Reliability flags, both established by Monte Carlo:
        %   * residual correlation makes LB anti-conservative as u/n grows;
        %   * the asymptotic chi2(h) needs n >> h (at n=100, h=20 the variance
        %     of Q was 52 against a nominal 40).
        ratio    = u / nc;
        reliable = (ratio <= 0.05) && (nc >= 20 * h);

        S.ljungBox = struct('Q', Q, 'p', p, 'h', h, 'reliable', reliable, ...
                            'uOverN', ratio);
        if ~reliable
            lsqvce.util.log.warn('lsqvce:systematic:lbUnreliable', ...
                ['Ljung-Box is anti-conservative here (u/n = %.3f, n = %d, ' ...
                 'h = %d). Least-squares residuals are correlated by ' ...
                 'construction -- cov(v) = s0^2*Qvv is not diagonal -- and the ' ...
                 'over-rejection grows with u/n. Treat its p-value as a screen, ' ...
                 'not a verdict, and lead with the binned test.'], ...
                 ratio, nc, h);
        end
    end
end

%% ----------------------------------------------------------- (c) binned means
names = fieldnames(o.covariates);
bins = struct('name', {}, 'z', {}, 'p', {}, 'worstBin', {}, 'worstZ', {}, ...
              'edges', {}, 'counts', {}, 'means', {});

for ii = 1:numel(names)
    cv = o.covariates.(names{ii});
    cv = cv(:);
    if numel(cv) ~= n
        error('lsqvce:systematic:badCovariate', ...
            ['Covariate ''%s'' has %d elements but there are %d ' ...
             'observations. It must be aligned with the observation order ' ...
             'the batch function produced.'], names{ii}, numel(cv), n);
    end
    cvc = cv(ctrl);
    good = isfinite(cvc);
    if sum(good) < 2 * o.nbins
        continue;
    end

    % Quantile edges: equal-count bins, so every bin has power. Equal-WIDTH
    % bins would leave the tails nearly empty and waste the Bonferroni budget
    % on bins that cannot reject anything.
    q  = linspace(0, 1, o.nbins + 1);
    ed = quantile(cvc(good), q);
    ed = unique(ed);
    nb = numel(ed) - 1;
    if nb < 2, continue; end
    ed(end) = ed(end) + eps(ed(end)) * 10;

    zb = zeros(nb, 1);  cnt = zeros(nb, 1);  mu = zeros(nb, 1);
    for b = 1:nb
        if b < nb
            m = good & cvc >= ed(b) & cvc < ed(b+1);
        else
            m = good & cvc >= ed(b) & cvc <= ed(b+1);
        end
        cnt(b) = sum(m);
        if cnt(b) < 2
            zb(b) = 0;  mu(b) = 0;
            continue;
        end
        mu(b) = mean(wc(m));
        % Under H0, w ~ N(0,1) so the bin mean has sd 1/sqrt(n_b).
        zb(b) = mu(b) * sqrt(cnt(b));
    end

    [~, ib] = max(abs(zb));
    praw = 2 * (1 - lsqvce.engine.stat.normCdf(abs(zb(ib))));
    padj = min(1, praw * nb);          % Bonferroni over the bins

    bins(end+1) = struct('name', names{ii}, 'z', zb, 'p', padj, ...
                         'worstBin', ib, 'worstZ', zb(ib), 'edges', ed, ...
                         'counts', cnt, 'means', mu); %#ok<AGROW>
end
S.bins = bins;

%% ------------------------------------------------------------------ verdict
ps = [];  lbl = {};
if isfinite(S.runs.p),     ps(end+1) = S.runs.p;     lbl{end+1} = 'runs'; end
if isfinite(S.ljungBox.p), ps(end+1) = S.ljungBox.p; lbl{end+1} = 'Ljung-Box'; end
for ii = 1:numel(bins)
    ps(end+1) = bins(ii).p; %#ok<AGROW>
    lbl{end+1} = sprintf('bins:%s', bins(ii).name); %#ok<AGROW>
end

if isempty(ps)
    S.detected = false;  S.worst = 'n/a';  S.pMin = NaN;
else
    [pmin, im] = min(ps);
    S.pMin     = pmin;
    S.worst    = lbl{im};
    S.detected = pmin < o.alpha;
end

%% ------------------------------------------------------------------- print
if ~o.verbose, return; end

lsqvce.util.log.msg('\n--- SYSTEMATIC ERROR TESTS  (alpha = %.3g)\n', o.alpha);

if isfinite(S.runs.p)
    lsqvce.util.log.msg('  runs test  : z = %+7.3f   p = %.3e   (%d runs, expected %.1f)%s\n', ...
        S.runs.z, S.runs.p, S.runs.nRuns, S.runs.expected, mark(S.runs.p, o.alpha));
end
if isfinite(S.ljungBox.p)
    rel = '';
    if ~S.ljungBox.reliable
        rel = '   [anti-conservative here: treat as a screen]';
    end
    lsqvce.util.log.msg('  Ljung-Box  : Q = %8.2f   p = %.3e   (h = %d)%s%s\n', ...
        S.ljungBox.Q, S.ljungBox.p, S.ljungBox.h, ...
        mark(S.ljungBox.p, o.alpha), rel);
end

if isempty(bins)
    lsqvce.util.log.msg(['  binned     : no covariates supplied.\n' ...
             '               This is the most powerful test and you have not\n' ...
             '               given it anything to work with. Pass every ordering\n' ...
             '               variable you have:\n' ...
             '                 ''covariates'', struct(''time'',t,''elev'',e,''azim'',az)\n']);
else
    for ii = 1:numel(bins)
        b = bins(ii);
        lsqvce.util.log.msg('  bins:%-8s: worst bin %d/%d, mean w = %+.4f (n=%d), z = %+.3f, p = %.3e%s\n', ...
            b.name, b.worstBin, numel(b.z), b.means(b.worstBin), ...
            b.counts(b.worstBin), b.worstZ, b.p, mark(b.p, o.alpha));
    end
end

lsqvce.util.log.msg('\n');
if S.detected
    lsqvce.util.log.msg(['  ==> SYSTEMATIC ERROR DETECTED (%s, p = %.2e).\n' ...
             '      The residuals carry structure that the model does not.\n' ...
             '      This is NOT an outlier problem and rejecting observations\n' ...
             '      will not fix it -- it will only hide it. Find the missing\n' ...
             '      parameter.\n'], S.worst, S.pMin);
    if ~isempty(bins)
        [~, iw] = min([bins.p]);
        b = bins(iw);
        lsqvce.util.log.msg(['      The signal is strongest against ''%s''. Plot w against\n' ...
                 '      it and look at the shape -- that shape IS the missing\n' ...
                 '      term.\n'], b.name);
    end
else
    lsqvce.util.log.msg(['  ==> No systematic structure detected at alpha = %.3g.\n' ...
             '      Note this is only as good as the covariates you supplied:\n' ...
             '      a systematic in a variable you did not pass is invisible\n' ...
             '      to every test here.\n'], o.alpha);
end
lsqvce.util.log.msg('\n');
end

function s = mark(p, a)
%MARK  Return a '<<< REJECT' tag when a p-value falls below the level a, else ''.
%   INPUT   p 1x1 double  a test p-value;  a 1x1 double  the significance level
%   OUTPUT  s 1xk char     the marker string ('' if not significant)
    if p < a, s = '   <<< REJECT'; else, s = ''; end
end
