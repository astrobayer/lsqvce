function [cfg, truth] = fixture_kbr_l1b(dataDir)
%FIXTURE_KBR_L1B  One synthetic day of GRACE-FO Level-1B files for KBR_VARIATIONAL.
%
%   [cfg, truth] = fixture_kbr_l1b(dataDir)
%
%   Writes GNV1B (C, D), KBR1B (Y), ACT1B and SCA1B (C, D) for MJD 60000 in
%   the GRACE-FO v04 ASCII layout of template/kbr_variational into dataDir,
%   from a TRUE field (the reference field + a planted correction on degrees
%   2..6) and a TRUE accelerometer calibration (bias, 1-3 % scale); returns
%   the cfg that processes them against the REFERENCE field and the a-priori
%   calibration (bias 0, scale 1), two 12 h arcs. The truth orbits come from
%   the template's own integrator (kit.propagate) -- the closed loop proves
%   the PIPELINE; the physics have their own oracles (J2, Kepler, re-
%   integration) in test_physOracle.
%
%   INPUT   dataDir  char  an existing, writable folder
%   OUTPUT  cfg      1x1 struct  the kbr_variational configuration
%           truth    1x1 struct  .dc (the planted correction, commonParams
%                    order), .accBias, .accScale (per satellite), .kit
%
%   ERRORS  none of its own (file errors of fopen pass through).
%
%   Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-122.
    day = 60000;  N = 6;  h = 5;
    gf = struct('GM', 3.986004415e14, 'R', 6378136.3, 'C', zeros(N + 1), 'S', zeros(N + 1));
    gf.C(1, 1) = 1;  gf.C(3, 1) = -4.841694573e-4;
    rs = RandStream('mt19937ar', 'Seed', 7);
    for n = 2:N
        for m = 0:n
            if ~(n == 2 && m == 0), gf.C(n + 1, m + 1) = 1e-5 / n^2 * randn(rs); end
            if m > 0, gf.S(n + 1, m + 1) = 1e-5 / n^2 * randn(rs); end
        end
    end
    cfg = struct('days', day, 'dataDir', dataDir, 'gravityField', gf, ...
                 'trf2crf', @(t) eye(3), 'degree', [2 N], 'arcsPerDay', 2, ...
                 'step', h, 'posSampling', 30);
    [~, kit] = kbr_variational(cfg);
    cfg.trf2crf = kit.eraRotation;

    % ---- the truth: field + planted correction, ACC calibration
    [nn, ~, ~] = kit.coeffIndex(N);
    est = find(nn >= 2);
    c = kit.field.c;
    dc = 2e-9 * randn(rs, numel(est), 1) ./ nn(est);
    cT = c;  cT(est) = cT(est) + dc;
    FT = kit.makeField(cT, N, [], gf.GM, gf.R);
    truth = struct('dc', dc, 'kit', kit);
    truth.accBias = struct('C', [3e-9; -2e-9; 1e-9], 'D', [-1e-9; 2.5e-9; 0.5e-9]);
    truth.accScale = struct('C', [1.02; 0.98; 1.01], 'D', [0.99; 1.03; 1.0]);

    % ---- the day: GPS seconds, a 5 s sample grid with one sample of margin
    t0 = (day - 51544.5) * 86400;
    tS = (t0 - h:h:t0 + 86400).';
    sep = 220e3 / (gf.R + 490e3);
    dateStr = char(datetime(day, 'ConvertFrom', 'modifiedjuliandate', 'Format', 'yyyy-MM-dd'));
    sat = 'CD';  states = cell(1, 2);
    for k = 1:2
        phase = 0.2 - (k - 1) * sep;
        [Q, acc] = attitudeAndAcc(tS - t0, phase);
        quat = zeros(numel(tS), 4);
        for i = 1:numel(tS), quat(i, :) = rot2quat(Q(:, :, i)); end
        accfun = kit.accModel(tS, acc, tS, quat, eye(3), 'srf2crf', [t0 t0 + 86400 - h], 10);
        y0 = leoState(gf, 0.4 - (k - 1) * sep);
        states{k} = kit.propagate(FT, kit.eraRotation, accfun, truth.accBias.(sat(k)), ...
                                  truth.accScale.(sat(k)), t0, y0, h, 86400 / h - 1);
        writeAcc(fullfile(dataDir, sprintf('ACT1B_%s_%s_04.txt', dateStr, sat(k))), tS, acc, sat(k));
        writeSca(fullfile(dataDir, sprintf('SCA1B_%s_%s_04.txt', dateStr, sat(k))), tS, quat, sat(k));
        writeGnv(fullfile(dataDir, sprintf('GNV1B_%s_%s_04.txt', dateStr, sat(k))), ...
                 t0 + (0:86400 / h - 1).' * h, states{k}, kit.eraRotation, sat(k));
    end
    dr = states{2}(:, 1:3) - states{1}(:, 1:3);  dv = states{2}(:, 4:6) - states{1}(:, 4:6);
    rr = sum(dr .* dv, 2) ./ sqrt(sum(dr.^2, 2));
    writeKbr(fullfile(dataDir, sprintf('KBR1B_%s_Y_04.txt', dateStr)), t0 + (0:86400 / h - 1).' * h, rr);
end

function [Q, acc] = attitudeAndAcc(t, phase)
%ATTITUDEANDACC  A slowly turning SRF -> CRF attitude and a drag-like ACC signal.
    w = 2 * pi / 5650;  n = numel(t);
    Q = zeros(3, 3, n);  acc = zeros(n, 3);
    T = [1 0 0; 0 cos(0.3) -sin(0.3); 0 sin(0.3) cos(0.3)];
    for i = 1:n
        th = w * t(i) + phase;
        Q(:, :, i) = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1] * T;
        acc(i, :) = [-1e-7 * (1 + 0.3 * cos(th)), 2e-8 * sin(th), 5e-8 * cos(2 * th)];
    end
end

function q = rot2quat(R)
%ROT2QUAT  Unit quaternion [q0 q1 q2 q3] (scalar first) of the rotation R.
    tr = trace(R);
    if tr > 0
        s = 2 * sqrt(1 + tr);
        q = [s / 4, (R(3, 2) - R(2, 3)) / s, (R(1, 3) - R(3, 1)) / s, (R(2, 1) - R(1, 2)) / s];
    else
        [~, k] = max(diag(R));
        switch k
            case 1
                s = 2 * sqrt(1 + R(1, 1) - R(2, 2) - R(3, 3));
                q = [(R(3, 2) - R(2, 3)) / s, s / 4, (R(1, 2) + R(2, 1)) / s, (R(1, 3) + R(3, 1)) / s];
            case 2
                s = 2 * sqrt(1 + R(2, 2) - R(1, 1) - R(3, 3));
                q = [(R(1, 3) - R(3, 1)) / s, (R(1, 2) + R(2, 1)) / s, s / 4, (R(2, 3) + R(3, 2)) / s];
            otherwise
                s = 2 * sqrt(1 + R(3, 3) - R(1, 1) - R(2, 2));
                q = [(R(2, 1) - R(1, 2)) / s, (R(1, 3) + R(3, 1)) / s, (R(2, 3) + R(3, 2)) / s, s / 4];
        end
    end
end

function y0 = leoState(gf, u)
%LEOSTATE  A circular 490 km, 89 deg orbit at argument of latitude u (CRF).
    a = gf.R + 490e3;  vc = sqrt(gf.GM / a);  i = deg2rad(89);
    y0 = [a * [cos(u); sin(u) * cos(i); sin(u) * sin(i)]; vc * [-sin(u); cos(u) * cos(i); cos(u) * sin(i)]];
end

function fid = openWithHeader(file, product)
    fid = fopen(file, 'w');
    fprintf(fid, 'header:\n  dimensions:\n    num_records: 0\n  global_attributes:\n    title: %s synthetic (lsqvce WP-122 test)\n# End of YAML header\n', product);
end

function writeGnv(file, t, states, rot, sat)
    fid = openWithHeader(file, 'GNV1B');  c = onCleanup(@() fclose(fid));
    for i = 1:numel(t)
        [R, Rd] = rot(t(i));
        r = R.' * states(i, 1:3).';  v = R.' * (states(i, 4:6).' - Rd * r);
        fprintf(fid, '%.1f %s E %.9f %.9f %.9f 0.01 0.01 0.01 %.12f %.12f %.12f 1e-5 1e-5 1e-5 00000000\n', ...
                t(i), sat, r, v);
    end
end

function writeKbr(file, t, rr)
    fid = openWithHeader(file, 'KBR1B');  c = onCleanup(@() fclose(fid));
    fprintf(fid, '%.1f 0 %.15e 0 0 0 0 0 0 0 0 1 1 1 1 00000000\n', [t, rr].');
end

function writeAcc(file, t, acc, sat)
    fid = openWithHeader(file, 'ACT1B');  c = onCleanup(@() fclose(fid));
    for i = 1:numel(t)
        fprintf(fid, '%.1f %s %.15e %.15e %.15e 0 0 0 0 0 0 00000000\n', t(i), sat, acc(i, :));
    end
end

function writeSca(file, t, quat, sat)
    fid = openWithHeader(file, 'SCA1B');  c = onCleanup(@() fclose(fid));
    for i = 1:numel(t)
        fprintf(fid, '%.1f %s 1 %.15f %.15f %.15f %.15f 0 00000000\n', t(i), sat, quat(i, :));
    end
end
