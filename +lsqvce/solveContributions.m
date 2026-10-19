function out = solveContributions(C, c, opt)
%LSQVCE.SOLVECONTRIBUTIONS  Solve a weighted combination of a run's days without re-reading them.
%
%   out = LSQVCE.SOLVECONTRIBUTIONS(C, c)
%   out = LSQVCE.SOLVECONTRIBUTIONS(C, c, 'Locals', true)
%
%   For a bootstrap over days (WP-158, issue #183). C = R.contributions()
%   holds every (system, day)'s reduced contribution of a fit; c gives each
%   day a multiplicity -- how often a replicate drew it (0, 1, 2, ...), or any
%   non-negative real for a Bayesian bootstrap. The combination's reduced
%   normal equations are the weighted sum of the days',
%       N* = sum_d c_d N_d,  n* = sum_d c_d n_d,  l'l* = sum_d c_d l'l_d,
%       dof* = sum_d c_d dof_d,
%   formed on the square root (each day's rows times sqrt(c_d)) on the QR
%   backend, and solved exactly as fit solves its accumulation: the common
%   prior added ONCE, the same rank rule and cancellation guards. A day drawn
%   twice contributes 2 N_d -- no relabelling around lsqvce:fit:duplicateDay.
%   One replicate costs O(D nc^2) for D days and nc common parameters,
%   independent of the number of observations per day.
%
%   THE WEIGHTS ARE FIXED. Every day carries the run's final weights: a noise
%   run's estimated variance components and a robust run's final (hard)
%   weights are NOT re-estimated per combination -- a bootstrap of the
%   estimator under fixed weights. out.weights says which weights those
%   were. A bootstrap that re-estimates them needs a full fit per replicate.
%   A robust run's sigma0 is divided by its truncation factor tau (C.truncation)
%   as fit's hard solution is: every combination of days is a sample of kept
%   epochs truncated alike.
%
%   INPUT
%     C   1x1 struct   R.contributions() of a linear fit
%     c   1 x D double multiplicity per day, aligned to C.days; finite, >= 0,
%                      at least one day > 0
%
%   OPTIONS (Name/Value)
%     'Locals'  logical  also back-substitute the local parameters of every
%                        day drawn (c_d > 0), with their covariance (default false)
%
%   OUTPUT
%     out.x       nc x 1   the common parameters
%     out.Cxx     nc x nc  their covariance (s0sq times the cofactor)
%     out.sigma0  1x1      sqrt(s0sq);  out.s0sq, out.dof, out.vtv as fit's solve
%     out.weights 1x1      C.weights: .estimate, .robust, .fixed = true
%     out.locals  struct array, with 'Locals' true, one per system and drawn
%                 day: .sysID, .day, .multiplicity, .localIdx (the system's
%                 local indices), .y, .Cyy; [] otherwise
%
%   EXAMPLE
%     R = lsqvce.fit(data, opts);
%     C = R.contributions();                       % one data pass
%     D = numel(C.days);  B = 100;  X = zeros(numel(R.x), B);
%     for b = 1:B
%         draw = randi(D, 1, D);                   % D days with replacement
%         c = accumarray(draw(:), 1, [D 1]).';
%         X(:, b) = lsqvce.solveContributions(C, c).x;
%     end
%     Cboot = cov(X.');                            % the bootstrap covariance of x
%
%   ACCURACY
%     With every multiplicity 1 the combination reproduces fit's x, Cxx and
%     sigma0 on the QR backend to 1e-8 relative, the issue's acceptance
%     bound (measured 2e-10 in x and Cxx, 2e-14 in sigma0 at cond(A) ~ 9e4:
%     the days' reduced rows are stacked and re-factorised where fit stacked
%     each system's); with a multiplicity vector it equals a full fit over the
%     expanded day list (a day drawn twice entered twice, under the same
%     fixed weights) to 1e-8 (measured 9e-10); the locals of a drawn day equal
%     fit's back-substitution to 1e-8 (measured 2e-10); a noise run and a
%     robust run (sigma0 with tau) reproduce their fit alike (test_fit
%     (testContributionsReproduceTheFit, testContributionsEqualAFitOverTheExpandedDays,
%     testContributionsReturnTheDrawnDaysLocals, testContributionsOfANoiseAndARobustRun)).
%     On the normal backend the agreement is that of the normal equations
%     themselves, cond(N) eps: at cond(N) ~ 7e9 fit's own normal solve is
%     5e-4 from the QR one (testContributionsOnTheNormalBackend uses a
%     well-conditioned fixture).
%
%   ERRORS
%     lsqvce:solveContributions:badContributions  C is not R.contributions()'s struct
%     lsqvce:solveContributions:badMultiplicity   c is not 1 x numel(C.days), finite,
%                                                 >= 0 with one day > 0
%     (and the solve's own: lsqvce:accum:priorNotPD, lsqvce:qr:*)
%
%   See also LSQVCE.FIT, LSQVCE.RESULTS
%
%   Created 2026-10-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-158 (issue #183).

arguments
    C (1,1) struct
    c double
    opt.Locals (1,1) logical = false
end

    need = {'nc', 'backend', 'prior', 'days', 'items', 'weights'};
    if ~all(isfield(C, need))
        error('lsqvce:solveContributions:badContributions', ...
              'C must be the struct R.contributions() returns (fields %s).', strjoin(need, ', '));
    end
    D = numel(C.days);
    if ~isequal(size(c), [1 D]) && ~isequal(size(c), [D 1])
        error('lsqvce:solveContributions:badMultiplicity', ...
              'c must hold one multiplicity per day of C.days (%d), got %d x %d.', D, size(c, 1), size(c, 2));
    end
    c = c(:).';
    if any(~isfinite(c)) || any(c < 0) || ~any(c > 0)
        error('lsqvce:solveContributions:badMultiplicity', ...
              'The multiplicities must be finite and non-negative, with at least one day drawn.');
    end
    acc = lsqvce.engine.BlockAccumulator.fromContributions(C, c);
    sol = acc.solve();
    tau = lsqvce.util.getfielddef(C, 'truncation', 1);   % a robust run's hard sigma0, as fit's
    sol.s0sq = sol.s0sq / tau;  sol.Cxx = sol.Cxx / tau;  sol.sigma0 = sqrt(sol.s0sq);
    out = struct('x', sol.x, 'Cxx', sol.Cxx, 'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'dof', sol.dof, ...
                 'vtv', lsqvce.util.getfielddef(sol, 'vtv', NaN), 'weights', C.weights, 'locals', []);
    if opt.Locals
        out.locals = drawnLocals(C, c, sol);
    end
end

function L = drawnLocals(C, c, sol)
%DRAWNLOCALS  The local parameters of every drawn day: y = R11 \ (r1 - R12 x),
%   Cyy = s0sq R11^-1 R11^-T + G Cxx G', G = R11 \ R12 (the normal backend:
%   Nll in place of R11'R11). Their covariance is the combination's.
    L = struct('sysID', {}, 'day', {}, 'multiplicity', {}, 'localIdx', {}, 'y', {}, 'Cyy', {});
    for k = 1:numel(C.items)
        it = C.items(k);
        w = c(C.days == it.day);
        if w == 0, continue, end
        if ~isempty(it.R11)
            G  = it.R11 \ it.R12;
            y  = it.R11 \ (it.r1 - it.R12 * sol.x);
            Ri = it.R11 \ eye(size(it.R11, 1));
            Cyy = sol.s0sq * (Ri * Ri.') + G * sol.Cxx * G.';
        elseif ~isempty(it.Nll)
            Ll = chol((it.Nll + it.Nll.') / 2, 'lower');
            G  = Ll.' \ (Ll \ it.Nlc);
            y  = Ll.' \ (Ll \ (it.nll - it.Nlc * sol.x));
            Li = Ll \ eye(size(Ll, 1));
            Cyy = sol.s0sq * (Li.' * Li) + G * sol.Cxx * G.';
        else
            G = zeros(0, numel(sol.x));  y = zeros(0, 1);  Cyy = zeros(0); %#ok<NASGU>
        end
        L(end+1) = struct('sysID', it.sysID, 'day', it.day, 'multiplicity', w, ...
                          'localIdx', it.localIdx, 'y', y, 'Cyy', (Cyy + Cyy.') / 2); %#ok<AGROW>
    end
end
