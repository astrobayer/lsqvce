function demo_combine()
%DEMO_COMBINE  Combining several satellite systems by elimination.
%
%   Three satellite systems observe a common gravity field. Each has its own
%   LOCAL parameters (arc states, accelerometer biases, ambiguities) that no
%   other system can see. Those locals are eliminated system by system, so the
%   running memory never exceeds O(uC^2 + max(uL)^2).
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

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
%     local parameters. That is fine -- they land in the same block.
%   * A parameter seen by two different SYSTEMS is COMMON, however local it
%     feels. Eliminating it as though it were local is not an approximation:
%     it is wrong, and the residuals will look fine.

uC  = 6;                       % common: gravity field coefficients
uL  = [12, 8, 15];             % locals: per-system arc states / biases
nb  = [40, 30, 50];            % batches per system
nk  = 25;                      % observations per batch

nSys = numel(uL);
fprintf('%d systems, %d common parameters, locals = %s\n', ...
        nSys, uC, mat2str(uL));
fprintf('total parameters = %d, of which %d are eliminated\n\n', ...
        uC + sum(uL), sum(uL));

%% ------------------------------------------------------------------------
% Build the data. In real work each batchFcn reads from disk / integrates an
% orbit; here we generate it.
% -------------------------------------------------------------------------
xC_true = randn(rs, uC, 1);
data = cell(nSys, 1);
for i = 1:nSys
    ni = nb(i) * nk;
    data{i}.xL_true = randn(rs, uL(i), 1);
    data{i}.AC = randn(rs, ni, uC);
    data{i}.AL = randn(rs, ni, uL(i));
    data{i}.sg = 0.5 + rand(rs, ni, 1);
    data{i}.p  = 1 ./ data{i}.sg.^2;
    data{i}.l  = data{i}.AC * xC_true ...
               + data{i}.AL * data{i}.xL_true ...
               + data{i}.sg .* randn(rs, ni, 1);
    data{i}.n  = ni;
end

% one batch source per system (a NESTED function: it sees `data` and `nk`)
    function [AC, AL, l, p] = srcOf(i, k)
        d = data{i};
        idx = (k-1)*nk + (1:nk);
        AC = d.AC(idx, :);
        AL = d.AL(idx, :);
        l  = d.l(idx);
        p  = d.p(idx);
    end

%% ------------------------------------------------------------------------
% COMBINE
% -------------------------------------------------------------------------
cmb = lsqvce.engine.Combine(uC, 'method', 'qr');

for i = 1:nSys
    cmb.addSystem(uL(i), @(k) srcOf(i, k), nb(i), ...
                  'name', sprintf('satellite %d', i));
end

sol = cmb.solve();
D   = cmb.diagnose(sol);

%% ------------------------------------------------------------------------
% THE INTEGRITY CHECK -- do this FIRST.
% -------------------------------------------------------------------------
% sum(r_i) = n - u_TOTAL, counting the ELIMINATED parameters.
%
% The locals were removed from the SOLVE, not from the MODEL: they still cost
% degrees of freedom. Two classic errors both violate this identity, and both
% lie in the OPTIMISTIC direction:
%
%   * using n - uC instead of n - u_total   -> s0^2 deflated by ~7%
%   * dropping the cross term in a'Qxx a    -> sum(r) inflated by ~7%,
%                                              every r_i too large, every MDB
%                                              too optimistic
fprintf('INTEGRITY\n');
fprintf('  n = %d,  u_total = %d  (uC = %d + uL = %d)\n', ...
        D.n, D.uTotal, D.uC, D.uL);
fprintf('  redundancy r = n - u_total = %d\n', D.r);
fprintf('  sum(r_i) = %.6f    rel err %.2e\n', D.sumR, D.sumRRelErr);
if D.sumRRelErr < 1e-8
    fprintf('  -> identity holds.\n\n');
else
    fprintf('  -> VIOLATED. Do not read on.\n\n');
end

%% ------------------------------------------------------------------------
% RESULTS
% -------------------------------------------------------------------------
fprintf('COMMON PARAMETERS\n');
fprintf('  %5s %12s %12s %12s %10s\n', 'par', 'estimate', 'truth', ...
        'error', 'sigma');
fprintf('  %s\n', repmat('-', 1, 56));
for j = 1:uC
    fprintf('  %5d %12.6f %12.6f %12.2e %10.4f\n', ...
            j, sol.xC(j), xC_true(j), sol.xC(j)-xC_true(j), ...
            D.s0 * sol.sigC(j));
end

fprintf('\n  s0 = %.4f   (a priori sigma0 = %.1f)   variance ratio = %.4f\n', ...
        D.s0, D.sigma0, D.varRatio);

%% ------------------------------------------------------------------------
% THE LOCALS COME BACK BY BACK-SUBSTITUTION
% -------------------------------------------------------------------------
% x_i = Nll^-1 (b_i - N_iC xC),  and their covariance is
%
%       Q_ii = Nll^-1 + S QCC S'
%
% The SECOND TERM IS NOT OPTIONAL. Nll^-1 alone is the covariance of the
% locals with the common parameters held FIXED -- but xC was ESTIMATED, and
% its uncertainty propagates into the locals. Drop it and the local sigmas
% come out too optimistic.
fprintf('\nLOCAL PARAMETERS (recovered by back-substitution)\n');
for i = 1:nSys
    e = norm(sol.xL{i} - data{i}.xL_true) / norm(data{i}.xL_true);
    fprintf('  satellite %d: %2d locals, ||x - truth||/||truth|| = %.2e, ', ...
            i, uL(i), e);
    fprintf('mean sigma = %.4f\n', D.s0 * mean(sol.sigL{i}));
end

%% ------------------------------------------------------------------------
% DIAGNOSTICS PER SYSTEM
% -------------------------------------------------------------------------
% The redundancy numbers use the FULL Qxx, including the cross terms between
% the common and the local blocks -- which is why they can be trusted.
fprintf('\nPER-SYSTEM DIAGNOSTICS\n');
fprintf('  %-12s %7s %9s %9s %9s\n', 'system', 'n', 'mean r_i', ...
        'min r_i', 'flagged');
fprintf('  %s\n', repmat('-', 1, 52));
for i = 1:nSys
    ri = D.rNum{i};  wi = D.w{i};
    nf = sum(abs(wi) > D.kCrit);
    fprintf('  %-12s %7d %9.4f %9.4f %9d\n', ...
            sprintf('satellite %d', i), numel(ri), mean(ri), min(ri), nf);
end
fprintf('\n  critical |w| = %.3f  (alpha_eff = %.1e, Bonferroni over n = %d)\n', ...
        D.kCrit, D.alphaEff, D.n);

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

function s = addComma(x)
    s = regexprep(sprintf('%d', x), '(\d)(?=(\d{3})+$)', '$1,');
end
