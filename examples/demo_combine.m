function demo_combine()
%DEMO_COMBINE  Combining several satellite systems by elimination.
%
%   Three satellite systems observe a common gravity field. Each has its own
%   LOCAL parameters (arc states, accelerometer biases, ambiguities) that no
%   other system can see. You DECLARE that structure -- one entry of
%   data.systems per system, its localParams stating the locals -- and
%   lsqvce.fit eliminates the locals system by system, so the running memory
%   never exceeds O(uC^2 + max(uL)^2).
%
%   INPUT   none
%   OUTPUT  none (prints the walk-through; ~seconds)
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- lsqvce.engine.Combine replaced by the
%   declarative multi-system data struct and one lsqvce.fit; verdict and
%   diagnostics from the Results. With assistance from Claude (Anthropic),
%   model Claude Fable 5.

clc;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------------
% THE DEFINITION OF "LOCAL" -- get this wrong and the answer is wrong.
% -------------------------------------------------------------------------
% A parameter is LOCAL if exactly ONE data set can see it.
% It is COMMON if two or more can.
%
% This is about WHICH DATA SETS OBSERVE IT, not about whether it feels like a
% nuisance parameter.
%
%   * Several SENSORS within one satellite system may all observe the same
%     local parameters. That is fine -- they land in the same block
%     (examples/demo_sensors.m shows the full hierarchy).
%   * A parameter seen by two different SYSTEMS is COMMON, however local it
%     feels. Eliminating it as though it were local is not an approximation:
%     it is wrong, and the residuals will look fine.

uC  = 6;                       % common: gravity field coefficients
uL  = [12, 8, 15];             % locals: per-system arc states / biases
nDay = [40, 30, 50];           % "days" (batches) per system
nk  = 25;                      % observations per day

nSys = numel(uL);
fprintf('%d systems, %d common parameters, locals = %s\n', ...
        nSys, uC, mat2str(uL));
fprintf('total parameters = %d, of which %d are eliminated\n\n', ...
        uC + sum(uL), sum(uL));

%% ------------------------------------------------------------------------
% Build the data. In real work each buildDay reads from disk / integrates an
% orbit; here we generate it once and let buildDay READ it (deterministic).
% -------------------------------------------------------------------------
xC_true = randn(rs, uC, 1);
gen = cell(nSys, 1);
for i = 1:nSys
    ni = nDay(i) * nk;
    gen{i}.xL_true = randn(rs, uL(i), 1);
    gen{i}.AC = randn(rs, ni, uC);
    gen{i}.AL = randn(rs, ni, uL(i));
    gen{i}.l  = gen{i}.AC * xC_true ...
              + gen{i}.AL * gen{i}.xL_true ...
              + randn(rs, ni, 1);            % unit noise: the lesson here is
                                             % STRUCTURE, not weighting
end

%% ------------------------------------------------------------------------
% DECLARE -- one entry of data.systems per satellite; localParams IS the
% statement "these unknowns are mine alone" (nl = rows).
% -------------------------------------------------------------------------
data = struct();
data.commonParams = num2cell((1:uC).');
for i = 1:nSys
    s = struct();
    s.buildDay    = makeBuildDay(gen{i}, nk);
    s.covariance.obs = struct('type', 'diag', 'sigma', 1);
    s.localParams = num2cell((1:uL(i)).');
    s.days        = 1:nDay(i);              % per-system schedule
    s.name        = sprintf('satellite %d', i);
    data.systems(i) = s;                    % id defaults to the index
end
R = lsqvce.fit(data, lsqvce.Options(Verbose=false));

%% ------------------------------------------------------------------------
% THE INTEGRITY CHECK -- do this FIRST.
% -------------------------------------------------------------------------
% The locals were removed from the SOLVE, not from the MODEL: they still cost
% degrees of freedom, so the redundancy is n - u_total, counting the
% ELIMINATED parameters. The classic errors (n - uC in the divisor; dropped
% cross terms in the redundancy numbers) both lie in the OPTIMISTIC
% direction; the report's degrees of freedom let you check.
n = sum(nDay) * nk;
fprintf('INTEGRITY\n');
fprintf('  n = %d,  u_total = %d  (uC = %d + sum uL = %d)\n', ...
        n, uC + sum(uL), uC, sum(uL));
fprintf('  redundancy r = n - u_total = %d,  reported dof = %d\n', ...
        n - uC - sum(uL), R.diag.r);

%% ------------------------------------------------------------------------
% RESULTS
% -------------------------------------------------------------------------
fprintf('\nCOMMON PARAMETERS\n');
fprintf('  %5s %12s %12s %12s %10s\n', 'par', 'estimate', 'truth', ...
        'error', 'sigma');
fprintf('  %s\n', repmat('-', 1, 56));
sig = sqrt(diag(R.sol.Cxx));
for j = 1:uC
    fprintf('  %5d %12.6f %12.6f %12.2e %10.4f\n', ...
            j, R.x(j), xC_true(j), R.x(j) - xC_true(j), sig(j));
end
fprintf('\n  sigma0 = %.4f   variance ratio = %.4f\n', R.sigma0, R.diag.varRatio);

R.report

%% ------------------------------------------------------------------------
% WHY BOTHER: the cost
% -------------------------------------------------------------------------
uTot = uC + sum(uL);
fprintf('\nCOST\n');
fprintf('  full normal matrix : %d x %d = %s entries\n', ...
        uTot, uTot, addComma(uTot^2));
fprintf('  reduced (common)   : %d x %d = %d entries\n', uC, uC, uC^2);
fprintf('  largest local block: %d x %d = %d\n', max(uL), max(uL), max(uL)^2);
fprintf('  -> peak memory is O(uC^2 + max(uL)^2) = %d, not O(u_total^2).\n', ...
        uC^2 + max(uL)^2);
fprintf('\n  Scale it up: 100 systems x 50 locals with 6 common parameters\n');
fprintf('     full     = %s entries\n', addComma((6 + 100*50)^2));
fprintf('     reduced  = %s entries    (%dx smaller)\n', ...
        addComma(6^2 + 50^2), round((6+100*50)^2 / (6^2 + 50^2)));

fprintf('\n');

end

%% ========================================================================
function f = makeBuildDay(g, nk)
%MAKEBUILDDAY  One system's buildDay over its pre-generated data (read-only,
%   hence deterministic). RAW rows, per the contract -- the toolbox whitens
%   from the declared spec.
    f = @buildDay;
    function D = buildDay(day, ~)
        rows = (day-1)*nk + (1:nk);
        D.A  = g.AC(rows, :);
        D.B  = g.AL(rows, :);
        D.l  = g.l(rows);
        D.t  = day*86400 + (0:nk-1).' * 5;   % seconds
    end
end

function s = addComma(x)
    s = regexprep(sprintf('%d', x), '(\d)(?=(\d{3})+$)', '$1,');
end
