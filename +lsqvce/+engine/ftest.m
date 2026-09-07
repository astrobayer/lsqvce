function F = ftest(solFull, solRed, opt)
%LSQVCE.ENGINE.FTEST  Is the extra parameter earning its place?
%
%   F = LSQVCE.ENGINE.FTEST(solFull, solRed) tests whether the additional parameters
%   in the FULL model give a statistically significant improvement in fit over
%   the REDUCED (nested) model.
%
%       F = [ (vPv_red - vPv_full) / (u_full - u_red) ] / s0sq_full
%           ~ F(u_full - u_red,  r_full)   under H0 (the extra params are zero)
%
%   The reduced model must be strictly NESTED inside the full one -- same
%   observations, same weights, a subset of the parameters. This is checked as
%   far as it can be: n must match, and vPv_red >= vPv_full (adding freedom can
%   only improve the fit, so a reduced model that fits BETTER means the two
%   solutions are not what you think they are).
%
%   WHY YOU NEED THIS
%   -----------------
%   The report tells you that a LOW variance factor may mean over-
%   parameterisation -- the model has enough freedom to absorb the noise, which
%   produces a beautiful fit and meaningless parameters. That warning is useless
%   without a way to act on it. This is the way to act on it.
%
%   The temptation is to look at the t-statistic x_j / sigma_j for each
%   parameter and drop the insignificant ones. That is fine for ONE parameter,
%   but it is wrong for several at once when they are correlated: two parameters
%   can each be individually insignificant while their JOINT contribution is
%   highly significant (and vice versa). The F-test asks the right question --
%   does this GROUP of parameters, taken together, earn its keep.
%
%   INPUT
%     solFull  1x1 struct  solution of the FULL model (from lsqvce.engine.adjust /
%                          Accumulator.solve); uses .vtv (or .vPv), .m/.n,
%                          .u/.rank and .sigma0/.s0 as available.
%     solRed   1x1 struct  solution of the nested REDUCED model, same fields.
%
%   NAME-VALUE
%     'alpha'    significance (default 0.05 -- model selection is not blunder
%                detection, and 0.001 here would retain almost nothing)
%     'verbose'  print (default true)
%
%   OUTPUT F
%     .F, .p, .df1, .df2, .crit
%     .significant   true if the extra parameters are justified
%     .dvPv          the improvement in vPv
%
%   EXAMPLE
%       % Is the extra parameter justified, or is it fitting noise?
%       Rfull = lsqvce.fit(struct('A', [A extra], 'l', l, 'p', p), Verbose=false);
%       Rred  = lsqvce.fit(struct('A', A,         'l', l, 'p', p), Verbose=false);
%
%       F = lsqvce.engine.ftest(Rfull.sol, Rred.sol);   % = Rfull.compare(Rred)
%
%       % F.significant   false for a pure-noise column (measured p = 0.548)
%       % F.p             the p-value
%
%       % The test assumes the REDUCED model is correct. If a systematic error
%       % is present, both models are wrong and the p-value means nothing --
%       % pass 'checkD' and it will refuse:
%       F = lsqvce.engine.ftest(sF, sR, 'checkD', S);
%
%   ERRORS
%     lsqvce:ftest:notNested       the two models are not nested, so the F ratio
%                                  has no meaning
%     lsqvce:ftest:differentData    they were not fitted to the same observations
%     lsqvce:ftest:noRedundancy     the fuller model has no degrees of freedom left
%     lsqvce:ftest:zeroVariance     the residual sum of squares is zero
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:ftest:invalidPremise                [warning]  ... ... ... ... ... 
%
%   See also LSQVCE.ENGINE.SYSTEMATIC, LSQVCE.ENGINE.ADJUST.
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
    solFull
    solRed
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 0.05
    opt.Verbose (1,1) logical = true
    opt.CheckD = []
end

o = struct('alpha', opt.Alpha, 'verbose', opt.Verbose, 'checkD', opt.CheckD);

uF = solFull.rank;   if isempty(uF), uF = solFull.u; end
uR = solRed.rank;    if isempty(uR), uR = solRed.u;  end

if solFull.n ~= solRed.n
    error('lsqvce:ftest:differentData', ...
        ['The two models were fitted to different numbers of observations ' ...
         '(%d and %d). An F-test compares two models on the SAME data; ' ...
         'otherwise the difference in vPv reflects the data, not the model.'], ...
         solFull.n, solRed.n);
end

df1 = uF - uR;
df2 = solFull.n - uF;

if df1 <= 0
    error('lsqvce:ftest:notNested', ...
        ['The "full" model has %d estimable parameters and the "reduced" one ' ...
         'has %d. The full model must have MORE.'], uF, uR);
end
if df2 <= 0
    error('lsqvce:ftest:noRedundancy', ...
        'The full model has no redundancy (n = %d, rank = %d).', solFull.n, uF);
end

dvPv = solRed.vPv - solFull.vPv;

if dvPv < 0
    if abs(dvPv) < 1e-8 * max(solFull.vPv, 1)
        dvPv = 0;      % roundoff
    else
        error('lsqvce:ftest:notNested', ...
            ['The REDUCED model fits BETTER than the full one ' ...
             '(vPv %.6e vs %.6e). That is impossible for genuinely nested ' ...
             'models -- adding parameters cannot worsen a least-squares fit.\n' ...
             'The two models are not nested, or they were not fitted to the ' ...
             'same data with the same weights, or one of them has not ' ...
             'converged.'], solRed.vPv, solFull.vPv);
    end
end

s0sqF = solFull.vPv / df2;
if s0sqF <= 0
    error('lsqvce:ftest:zeroVariance', 'The full model has zero residual.');
end

Fstat = (dvPv / df1) / s0sqF;
crit  = lsqvce.engine.stat.fInv(1 - o.alpha, df1, df2);
% p-value from the F CDF, via the incomplete beta
xb = (df1 * Fstat) / (df1 * Fstat + df2);
pval = 1 - betainc(xb, df1/2, df2/2);

F = struct('F', Fstat, 'p', pval, 'df1', df1, 'df2', df2, ...
           'crit', crit, 'alpha', o.alpha, 'dvPv', dvPv, ...
           'significant', Fstat > crit, ...
           'vPvFull', solFull.vPv, 'vPvRed', solRed.vPv);

% ---------------------------------------------------------------------------
% VALIDITY GUARD
% ---------------------------------------------------------------------------
% The F-test asks "does the extra parameter improve the fit MORE THAN CHANCE?"
% and answers it by comparing dvPv against the full model's residual variance.
% That comparison is only meaningful if the residual variance is NOISE.
%
% If a systematic error is present, or the data are contaminated, then s0sqF
% measures the systematic error rather than the noise, and the F-statistic is
% the ratio of two things that are both wrong. It will happily declare a
% parameter "significant" because it soaks up part of an unmodelled signal --
% which is precisely how a spurious parameter gets into a model and stays there.
F.valid = true;
F.validityWarning = '';
if ~isempty(o.checkD) && isstruct(o.checkD)
    bad = '';
    if isfield(o.checkD, 'contamLevel') && ...
       any(strcmp(o.checkD.contamLevel, {'moderate','severe','refuse','breakdown'}))
        bad = sprintf('contamination (level "%s")', o.checkD.contamLevel);
    elseif isfield(o.checkD, 'detected') && o.checkD.detected
        bad = 'an unmodelled systematic error';
    end
    if ~isempty(bad)
        F.valid = false;
        F.validityWarning = bad;
        lsqvce.util.log.warn('lsqvce:ftest:invalidPremise', ...
            ['ftest: the residuals show %s, so this F-test is NOT valid.\n' ...
             'The F-statistic divides the improvement in fit by the full ' ...
             'model''s residual variance, which is only meaningful when that ' ...
             'variance is NOISE. Here it is measuring the model error itself, ' ...
             'so the test compares two wrong things.\n' ...
             'A parameter can come out "significant" simply because it soaks up ' ...
             'part of an unmodelled signal -- which is exactly how a spurious ' ...
             'parameter enters a model and stays there. Fix the model first.'], ...
             bad);
    end
end

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\n--- F-TEST: are the %d extra parameter(s) justified?\n', df1);
    lsqvce.util.log.msg('  vPv (reduced, %d par) = %.6e\n', uR, solRed.vPv);
    lsqvce.util.log.msg('  vPv (full,    %d par) = %.6e\n', uF, solFull.vPv);
    lsqvce.util.log.msg('  improvement          = %.6e\n', dvPv);
    lsqvce.util.log.msg('  F = %.4f   ~ F(%d, %d)   critical = %.4f   p = %.4e\n', ...
            Fstat, df1, df2, crit, pval);
    lsqvce.util.log.msg('\n');
    if F.significant
        lsqvce.util.log.msg(['  ==> SIGNIFICANT. The extra parameter(s) reduce the residual\n' ...
                 '      by more than chance would. Keep them.\n']);
    else
        lsqvce.util.log.msg(['  ==> NOT significant. The extra parameter(s) buy an\n' ...
                 '      improvement no bigger than noise. They are absorbing\n' ...
                 '      random error, not signal -- which makes the fit look\n' ...
                 '      better while making the parameters MEAN less. Drop them.\n']);
    end
    lsqvce.util.log.msg('\n');
end
end
