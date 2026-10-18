function T = wp156_selection_measure(outFile, R)
%WP156_SELECTION_MEASURE  WP-156's controlled study: which statistic a
%   delete-then-divide-by-tau robust variance estimator must reject on.
%
%   T = WP156_SELECTION_MEASURE(outFile, R)
%
%   A MEASUREMENT, not a test: it asserts nothing and prints; HANDOVER
%   12.91 quotes its output. The noise law is KNOWN, so only the selection
%   is under test. Per rep 40 independent arcs of 256 epochs, AR(1) phi,
%   marginal sigma 1, a local bias + trend per arc (GLS under the true law),
%   optionally blunders (fraction eps, size b marginal sigma, random sign).
%   Epochs are rejected beyond kRej centred-MAD units of a statistic and
%   deleted exactly (rows removed); level = Omega / r / tau(kRej), tau the
%   white Gaussian truncation factor (adjustVCEjoint>truncationFactor).
%   The statistics:
%     raw  the raw residual v_i -- what robustSweep rejects on today
%          (rejection simultaneous, iterated to a fixed point);
%     w    Baarda's w for correlated observations, (P v)_i / sqrt((P Qvv
%          P)_ii) on kept epochs, the predictive statistic on deleted ones
%          (simultaneous, to a fixed point);
%     ws   w with iterative data snooping per arc: the largest |w| / s >
%          kRej deleted, recomputed, until none is left; s from the
%          undeleted system;
%     ws2  ws with s re-estimated from the final system (snooped epochs
%          +-Inf) and the snooping redone until s settles.
%   Every row also reports option 2 of review 12.84 N3 on the same run:
%   the level divided by an analytic coloured tau, (1 - mean_i[rho_i^2
%   E(x^2; |x| > k) + (1 - rho_i^2) p]) / (1 - p), rho_i^2 = c_i / sigma_ii
%   (meaningful on the raw rows only).
%   Why w: deleting epoch i removes exactly w_i^2 from Omega, so a rule that
%   selects on w truncates the quantity it removes and tau is the white
%   value for ANY law; the raw residual is correlated with w_i only by
%   rho_i = sqrt(c_i / sigma_ii), and the level reads high by the missing
%   truncation (review 12.84 N3).
%
%   INPUT   outFile  char  the .mat file the tables are saved to after
%                          every configuration ('' for none)
%           R        1x1   reps per clean configuration (default 400; the
%                          contaminated ones run R * 3/8)
%   OUTPUT  T  struct array, one per configuration: .name, per rep .lev
%              (level / truth), .levM (option 2), .nRej, .det (blunders
%              rejected), .fa (clean epochs rejected), .xe / .xo (squared
%              error of the local parameters, run / oracle), .iters
%
%   See also wp115_measure
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-156, the design study.
    if nargin < 1, outFile = ''; end
    if nargin < 2 || isempty(R), R = 400; end
    cfgs = struct('phi', {}, 'eps', {}, 'b', {}, 'kRej', {}, 'stat', {}, 'name', {});
    for ph = [0 0.5 0.9 0.98]
        for st = {'raw', 'w', 'ws', 'ws2'}
            cfgs(end+1) = struct('phi', ph, 'eps', 0, 'b', 0, 'kRej', 3, 'stat', st{1}, ...
                'name', sprintf('clean_%s_phi%.2f', st{1}, ph)); %#ok<AGROW>
        end
    end
    T = runSet(outFile, cfgs, R);
    cfgs = cfgs([]);
    for ph = [0 0.9 0.98]
        for b = [2 4 8]
            for st = {'raw', 'w', 'ws', 'ws2'}
                cfgs(end+1) = struct('phi', ph, 'eps', 0.02, 'b', b, 'kRej', 3, 'stat', st{1}, ...
                    'name', sprintf('cont_%s_phi%.2f_b%d', st{1}, ph, b)); %#ok<AGROW>
            end
        end
    end
    T = [T, runSet(outFile, cfgs, max(1, round(R * 3 / 8)))];
end

function T = runSet(outFile, cfgs, R)
    T = struct('name', {}, 'lev', {}, 'levM', {}, 'nRej', {}, 'det', {}, 'fa', {}, 'xe', {}, 'xo', {}, 'iters', {});
    for ci = 1:numel(cfgs)
        c = cfgs(ci);
        lev = nan(R, 1);  levM = lev;  nRej = lev;  det = lev;  fa = lev;  xe = lev;  xo = lev;  it = lev;
        for rr = 1:R
            [lev(rr), levM(rr), nRej(rr), det(rr), fa(rr), xe(rr), xo(rr), it(rr)] = oneRep(c, rr);
        end
        T(ci) = struct('name', c.name, 'lev', lev, 'levM', levM, 'nRej', nRej, 'det', det, 'fa', fa, ...
                       'xe', xe, 'xo', xo, 'iters', it);
        if ~isempty(outFile), save(outFile, 'T', 'cfgs'); end
        fprintf('%-24s level %.4f +- %.4f (option 2: %.4f)  rejected %.4f  blunders caught %.3f  clean rejected %.5f  x eff %.3f\n', ...
            c.name, mean(lev), std(lev) / sqrt(R), mean(levM), mean(nRej), mean(det, 'omitnan'), mean(fa), ...
            mean(xe) / mean(xo));
    end
end

function [lev, levM, nRej, det, fa, xe, xo, iters] = oneRep(c, rr)
    nArc = 40;  m = 256;
    rs = RandStream('mt19937ar', 'Seed', 1000 * rr + round(1000 * c.phi) + 7);
    S = toeplitz(c.phi .^ (0:m - 1));
    L = chol(S, 'lower');
    t = ((1:m).' - (m + 1) / 2) / m;
    A = [ones(m, 1), t];
    xt = [1; 2];
    l = cell(nArc, 1);  bl = cell(nArc, 1);  lc = cell(nArc, 1);
    for a = 1:nArc
        e = L * randn(rs, m, 1);
        lc{a} = A * xt + e;
        bl{a} = rand(rs, m, 1) < c.eps;
        sg = sign(randn(rs, m, 1));
        l{a} = lc{a} + c.b * sg .* bl{a};
    end
    % oracle: clean data, all epochs, true law
    Pi = inv(S);  Pi = (Pi + Pi.') / 2;
    xo = 0;
    for a = 1:nArc
        x = (A.' * Pi * A) \ (A.' * Pi * lc{a});
        xo = xo + sum((x - xt).^2);
    end
    D = cell(nArc, 1);  for a = 1:nArc, D{a} = false(m, 1); end
    iters = 0;
    if any(strcmp(c.stat, {'ws', 'ws2'}))
        [D, om, r, xe, iters] = snoop(S, A, l, xt, c.kRej, nArc, m, strcmp(c.stat, 'ws2'));
        pass = 0;
    end
    for pass = 1:15 * ~any(strcmp(c.stat, {'ws', 'ws2'}))
        st = cell(nArc, 1);  om = 0;  r = 0;  xe = 0;
        for a = 1:nArc
            [st{a}, oa, ra, x] = arcStat(S, A, l{a}, D{a}, c.stat);
            om = om + oa;  r = r + ra;  xe = xe + sum((x - xt).^2);
        end
        u = vertcat(st{:});
        s = 1.4826 * median(abs(u - median(u)));
        rejAll = abs(u / s) > c.kRej;
        Dn = mat2cell(rejAll, m * ones(nArc, 1), 1);
        iters = pass;
        if isequal(Dn, D), break, end
        D = Dn;
    end
    % the fixed point's level (D is what this pass's om / r used)
    k = c.kRej;  p = erfc(k / sqrt(2));
    tau = 1 - 2 * k * normpdfL(k) / (1 - p);
    lev = om / r / tau;
    % option 2's analytic tau (from the conditional-variance share of each epoch,
    % no deletion): mean over epochs
    ci = 1 ./ diag(Pi);  rho2 = ci ./ diag(S);
    m2 = 2 * (k * normpdfL(k)) + p;                   % E[x^2; |x| > k]
    tauM = (1 - mean(rho2 * m2 + (1 - rho2) * p)) / (1 - p);
    levM = om / r / tauM;
    Dall = vertcat(D{:});  Ball = vertcat(bl{:});
    nRej = mean(Dall);
    det = NaN;  if any(Ball), det = mean(Dall(Ball)); end
    fa = mean(Dall(~Ball));
end

function [st, om, r, x] = arcStat(S, A, l, D, kind)
    K = ~D;
    SK = S(K, K);  W = chol(SK);                    % upper, SK = W' W
    Ab = W.' \ A(K, :);  lb = W.' \ l(K);
    [Q, Rr] = qr(Ab, 0);
    x = Rr \ (Q.' * lb);
    vb = lb - Ab * x;
    om = vb.' * vb;  r = nnz(K) - size(A, 2);
    v = l - A * x;
    m = numel(l);  st = zeros(m, 1);
    switch kind
        case 'raw'
            st = v;
        case 'w'
            Wi = inv(W);                            % SK^-1 = Wi Wi'
            Pv = Wi * vb;                           % SK^-1 v_K
            G = Wi * (Ab / Rr);                     % SK^-1 A_K R^-1 (Qxx = R^-1 R^-T)
            dP = sum(Wi.^2, 2) - sum(G.^2, 2);
            st(K) = Pv ./ sqrt(dP);
            % deleted epochs: the predictive statistic (what deleting them removes)
            idx = find(D);
            if ~isempty(idx)
                Z = W.' \ S(K, idx);                % W^-T S_Ki
                cvar = 1 - sum(Z.^2, 1).';          % conditional variances (sigma_ii = 1)
                ui = v(idx) - Z.' * (W.' \ v(K));
                Hi = A(idx, :) - Z.' * Ab;          % the prediction's design part
                vx = sum((Hi / Rr).^2, 2);          % its variance through x
                st(idx) = ui ./ sqrt(cvar + vx);
            end
    end
end

function y = normpdfL(x)
    y = exp(-x.^2 / 2) / sqrt(2 * pi);
end

function [D, om, r, xe, nDel] = snoop(S, A, l, xt, k, nArc, m, rescale)
%SNOOP  Baarda's iterative data snooping per arc (the arcs carry only local
%   parameters, so each arc snoops on its own): w of every epoch with nothing
%   deleted, the stream's centred MAD from those; then per arc the largest
%   |w| / s > k deleted, recomputed, until none is left. RESCALE: the MAD
%   again from the final system (kept w, snooped epochs +-Inf), and the
%   snooping redone from scratch with it, until the scale settles.
    w0 = cell(nArc, 1);
    for a = 1:nArc, w0{a} = arcStat(S, A, l{a}, false(m, 1), 'w'); end
    u = vertcat(w0{:});  s = 1.4826 * median(abs(u - median(u)));
    for round = 1:(1 + 9 * rescale)
        D = cell(nArc, 1);  om = 0;  r = 0;  xe = 0;  nDel = 0;  wf = cell(nArc, 1);
        for a = 1:nArc
            Da = false(m, 1);  w = w0{a};
            while true
                ww = abs(w);  ww(Da) = 0;
                [mx, j] = max(ww);
                if mx <= k * s, break, end
                Da(j) = true;  nDel = nDel + 1;
                w = arcStat(S, A, l{a}, Da, 'w');
            end
            [w, oa, ra, x] = arcStat(S, A, l{a}, Da, 'w');
            w(Da) = Inf;  wf{a} = w;
            om = om + oa;  r = r + ra;  xe = xe + sum((x - xt).^2);  D{a} = Da;
        end
        if ~rescale, break, end
        u = vertcat(wf{:});  sNew = 1.4826 * median(abs(u - median(u)));
        if abs(sNew / s - 1) < 1e-9, break, end
        s = sNew;
    end
end
