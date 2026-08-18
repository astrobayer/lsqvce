function bench_qr_vs_normal(varargin)
%BENCH_QR_VS_NORMAL  Measure the QR / normal-equation trade-off on YOUR machine.
%
%   The numbers quoted in the class documentation were measured with NumPy on
%   a container. MATLAB uses the same BLAS/LAPACK, so the RATIO should carry
%   over -- but the absolute times will not, and your BLAS threading, your u,
%   and your batch sizes are not mine. Measure.
%
%   BENCH_QR_VS_NORMAL()              a default sweep
%   BENCH_QR_VS_NORMAL('u', [10 50 100], 'n', 2e5, 'nk', 5000)
%
%   WHAT YOU SHOULD SEE
%     * QR costs roughly 2-9x the time, growing with u, not with n.
%     * Both paths are LINEAR in n.
%     * The QR cost is FLAT in the conditioning -- you do not pay more for QR
%       precisely when you need it more.
%
%   THE DECISION IS NOT ABOUT SPEED, IT IS ABOUT cond(N)
%     cond(N) < 1e8      normal loses ~8 of 16 digits. Fine. Take the speedup.
%     cond(N) 1e8..1e12  borderline. Use QR.
%     cond(N) > 1e12     normal has ~4 digits left. QR is not optional.
%
%   Since cond(N) = cond(A)^2, one cheap 'normal' run tells you which regime
%   you are in. That is the recipe: solve once with 'normal', read cond(N) from
%   the report, and re-run with 'qr' only if you must.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

ip = inputParser;
addParameter(ip, 'u',    [10 30 50 100], @isnumeric);
addParameter(ip, 'n',    100000,          @isscalar);
addParameter(ip, 'nk',   5000,            @isscalar);
addParameter(ip, 'reps', 3,               @isscalar);
parse(ip, varargin{:});
o = ip.Results;

rs = RandStream('mt19937ar', 'Seed', 1);

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  QR vs NORMAL EQUATIONS  (n = %d, batch = %d)\n', o.n, o.nk);
fprintf('%s\n', repmat('=', 1, 70));
fprintf('  %5s %11s %11s %10s %12s\n', ...
        'u', 'normal (s)', 'QR (s)', 'QR/normal', 'agreement');
fprintf('  %s\n', repmat('-', 1, 56));

for u = o.u(:).'
    A = randn(rs, o.n, u);
    l = randn(rs, o.n, 1);
    p = 0.5 + rand(rs, o.n, 1);

    e = [0:o.nk:o.n, o.n];
    e = unique(e);
    nb = numel(e) - 1;
    src = @(k) deal(A(e(k)+1:e(k+1), :), l(e(k)+1:e(k+1)), p(e(k)+1:e(k+1)));

    tN = inf;  tQ = inf;
    for rep = 1:o.reps
        t = tic;
        accN = lsqvce.Accumulator(u, 'method', 'normal');
        for k = 1:nb, [Ak,lk,pk] = src(k); accN.add(Ak,lk,pk); end
        solN = accN.solve();
        tN = min(tN, toc(t));

        t = tic;
        accQ = lsqvce.Accumulator(u, 'method', 'qr');
        for k = 1:nb, [Ak,lk,pk] = src(k); accQ.add(Ak,lk,pk); end
        solQ = accQ.solve();
        tQ = min(tQ, toc(t));
    end

    agree = norm(solN.x - solQ.x) / max(norm(solN.x), realmin);
    fprintf('  %5d %11.4f %11.4f %9.2fx %12.1e\n', u, tN, tQ, tQ/tN, agree);
end

fprintf('\n');
fprintf('  The two paths AGREE to ~1e-14 on a well-conditioned problem.\n');
fprintf('  They diverge only when cond(N) is large -- and that is exactly\n');
fprintf('  when the normal-equation answer is the wrong one.\n\n');

%% -------------------------------------------------------------------------
% The part that actually decides it: accuracy vs conditioning.
% -------------------------------------------------------------------------
fprintf('%s\n', repmat('=', 1, 70));
fprintf('  WHEN IS THE ACCURACY WORTH THE TIME?\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('  A near-collinear column pair, swept. Reference: SVD (does not\n');
fprintf('  form N, so it is not subject to the squaring).\n\n');
fprintf('  %8s %10s %12s %12s %s\n', ...
        'eps', 'cond(N)', 'normal err', 'QR err', 'verdict');
fprintf('  %s\n', repmat('-', 1, 66));

u  = 10;
n  = 20000;
nk = 5000;
for eps = [1e-1 1e-3 1e-5 1e-7 1e-8]
    A = randn(rs, n, u);
    A(:,2) = A(:,1) + eps*randn(rs, n, 1);      % near-collinear
    xt = randn(rs, u, 1);
    l  = A*xt + 1e-3*randn(rs, n, 1);
    p  = ones(n,1);

    % reference: SVD on the whitened system
    [U,Sv,Vt] = svd(sqrt(p).*A, 'econ');
    xr = Vt * ((U' * (sqrt(p).*l)) ./ diag(Sv));

    e = unique([0:nk:n, n]);
    nb = numel(e) - 1;
    src = @(k) deal(A(e(k)+1:e(k+1), :), l(e(k)+1:e(k+1)), p(e(k)+1:e(k+1)));

    accN = lsqvce.Accumulator(u, 'method', 'normal');
    accQ = lsqvce.Accumulator(u, 'method', 'qr');
    for k = 1:nb
        [Ak,lk,pk] = src(k);
        accN.add(Ak,lk,pk);  accQ.add(Ak,lk,pk);
    end
    W = warning('off', 'all');
    sN = accN.solve();  sQ = accQ.solve();
    warning(W);

    eN = norm(sN.x - xr)/norm(xr);
    eQ = norm(sQ.x - xr)/norm(xr);

    if eN < 1e-9
        v = 'normal is fine';
    elseif eN < 1e-6
        v = 'normal ok, QR safer';
    elseif eN < 1e-3
        v = 'QR RECOMMENDED';
    else
        v = 'QR MANDATORY';
    end
    fprintf('  %8.0e %10.1e %12.2e %12.2e %s\n', ...
            eps, sQ.condN, eN, eQ, v);
end

fprintf(['\n  READ cond(N) FROM YOUR OWN REPORT. That single number tells you\n' ...
         '  which column you are in -- and therefore whether the 2-9x cost of\n' ...
         '  QR is buying you anything at all.\n\n']);
end
