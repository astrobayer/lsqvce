function rec = timePair(fA, fB, opt)
%TIMEPAIR  Time two callables against each other and return a ratio record.
%   rec = timePair(fA, fB)
%   rec = timePair(fA, fB, 'Name', value, ...)
%
%   THE one timing instrument of the suite (R5 M0, WP-78; VALIDATION_GUIDE
%   Part 4). Every speed budget is a RATIO between two sides doing the same
%   work on the same arrays, and every ratio is measured here and nowhere
%   else, so the next repair to the method is applied once. What it does,
%   and what no test does by hand:
%     1. one untimed warm-up call per side;
%     2. BOTH sides are timed inside ONE repeat -- never A to completion,
%        then B (drift lands on the later one; measured on the source
%        project: 8.3 inside the full runner vs 4.8-5.5 alone);
%     3. the order of the two sides ROTATES between repeats (a fixed order
%        leaves a small, always-identically-signed bias);
%     4. Repeats = 15 by default (measured spread over 10 trials:
%        5 -> 31 %, 9 -> 14 %, 15 -> 10 %, 21 -> 4.4 %);
%     5. the ratio is formed PER REPEAT, then the MEDIAN of those is the
%        statistic -- never a ratio of medians;
%     6. the record carries the per-repeat min..max band, both absolutes and
%        a machine tag, so a red budget can be told from a bad minute.
%   A budget is then a one-line assertion against rec.ratio with the
%   provenance of its expected value in a tol: comment at the site.
%
%   INPUT
%     fA, fB   1x1 function_handle  zero-argument callables (wrap arguments
%              in the handle). Side A is the NUMERATOR: rec.ratio = tA / tB.
%   NAME-VALUE
%     'Repeats'  1x1 double  >= 3, default 15   timed repeats
%     'Inner'    1x1 double  >= 1, default 1    calls per timing of one side
%                                              (raise it below ~1 ms, where
%                                              the timer is what is measured)
%     'Label'    1xK char/string, default ''   carried into the record
%   OUTPUT
%     rec  1x1 struct
%       .ratio    1x1 double   median over repeats of (tA/tB)
%       .band     1x2 double   [min max] of the per-repeat ratios
%       .tA, .tB  1x1 double   median seconds per call of each side
%       .ratios   Rx1 double   the per-repeat ratios (R = Repeats)
%       .repeats, .inner       the settings used
%       .label    char         the label
%       .machine  char         computer() + MATLAB release, for the log
%
%   ERRORS
%     lsqvce:timePair:badHandle   fA or fB is not a function handle
%     MATLAB:validators:*         Repeats / Inner outside their ranges
%
%   ACCURACY
%     Timing, not numerics: the record is EVIDENCE with a band, never a
%     number to compare bitwise. Two identical callables give a ratio near
%     1 with a band that reflects the machine's jitter (test_instruments
%     pins ratio in [0.5, 2] and band containing 1 on the identity pair).
%
%   EXAMPLE
%     rec = timePair(@() lsqvce.fit(data, opts), ...
%                    @() lsqvce.engine.adjustBlock(L, days, 'Verbose', false), ...
%                    'Label', 'fit dispatch');
%     % tol: expected 1.00 (dispatch adds validation only), budget 3x -> 1.25
%     % measured 2026-08-31 L4 postflight; slow is the risk => <=
%     verifyLessThanOrEqual(tc, rec.ratio, 1.25, ...
%         sprintf('%s: ratio %.2f band [%.2f %.2f]', rec.label, rec.ratio, rec.band));
%
%   See also TIMEIT, RUN_ALL_LSQVCE, GATEVERDICT
%
%   Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0, WP-78.

arguments
    fA
    fB
    opt.Repeats (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opt.Repeats, 3)} = 15
    opt.Inner   (1,1) double {mustBeInteger, mustBePositive} = 1
    opt.Label   (1,:) {mustBeTextScalar} = ''
end
    if ~isa(fA, 'function_handle') || ~isa(fB, 'function_handle')
        error('lsqvce:timePair:badHandle', ...
              'fA and fB must be zero-argument function handles.');
    end
    R = opt.Repeats;  n = opt.Inner;

    fA();  fB();                                   % 1. warm-up, untimed

    tA = zeros(R, 1);  tB = zeros(R, 1);
    for r = 1:R
        if mod(r, 2) == 1                          % 3. rotation
            tA(r) = timeOne(fA, n);  tB(r) = timeOne(fB, n);
        else
            tB(r) = timeOne(fB, n);  tA(r) = timeOne(fA, n);
        end
    end
    ratios = tA ./ tB;                             % 5. per repeat

    rec = struct( ...
        'ratio',   median(ratios), ...
        'band',    [min(ratios) max(ratios)], ...
        'tA',      median(tA), ...
        'tB',      median(tB), ...
        'ratios',  ratios, ...
        'repeats', R, ...
        'inner',   n, ...
        'label',   char(opt.Label), ...
        'machine', sprintf('%s %s', computer(), version('-release')));
end

function t = timeOne(f, n)
%TIMEONE  Seconds per call of f over n back-to-back calls (one timed block).
    t0 = tic;
    for k = 1:n
        f();
    end
    t = toc(t0) / n;
end
