function R = demo_gravity(opt)
%DEMO_GRAVITY  GRACE-like gravity recovery on a PHYSICALLY STRUCTURED design.
%
%   R = DEMO_GRAVITY()
%   R = DEMO_GRAVITY('MaxOuter', 6, 'Verbose', true)
%   R = DEMO_GRAVITY('Plots', true, 'ExportDir', '/tmp/figs')
%
%   The other demos exercise the algebra on random designs. This one runs the
%   same loop on a design whose every column is a real spherical-harmonic
%   potential partial evaluated along a sampled near-polar orbit, so the
%   normal matrix carries what a production gravity recovery carries: order-
%   block correlation, polar-gap conditioning, degree attenuation. Nothing in
%   the design is random -- only the noise realisation is, and that is seeded.
%
%   WHAT IT SHOWS, in order
%     1. THE PROBLEM. Two satellites (GRACE-C and a trailing GRACE-D) and the
%        inter-satellite potential difference read as an energy-balance KBR
%        stream, over five days of a RAGGED arc schedule: a different number
%        of arcs per day, missing arcs by ID, and a different length for every
%        arc. Each system carries its OWN number of per-arc empirical
%        parameters (2, 3 and 1), so the block-arrow solver meets three
%        different local block sizes in one adjustment.
%     2. THE WEIGHTS. The covariance declared to the toolbox has the right AR
%        shape but unit sigma on all three streams, while the observations
%        were generated with sigmas 1.0 / 2.0 / 0.5. VCE has to find that,
%        under real column correlation rather than the benign Gaussian
%        conditions of the other demos.
%     3. THE FIELD. The coefficients were planted by a closed formula plus ONE
%        anomaly at degree 5, order 3. The demo reports how far the estimate
%        sits from the truth IN UNITS OF ITS OWN FORMAL ERROR -- the only
%        honest way to read a gravity solution -- and the error degree
%        variances against their formal expectation.
%
%   OPTIONS
%     'MaxOuter'   1x1 double   outer VCE iterations           (default 4)
%     'Verbose'    1x1 logical  print the toolbox's own report (default false)
%     'Plots'      1x1 logical  draw the two summary figures   (default false)
%     'ExportDir'  1xK char     if non-empty, save the figures there as PNG
%                                                              (default '')
%
%   OUTPUT
%     R   1x1 lsqvce.Results   the full result: R.x, R.sigma0, R.sol.Cxx,
%                              R.noise (the recovered stream models), R.meta
%                              (run provenance), R.report() for the verdict.
%
%   ERRORS
%     lsqvce:demoGravity:noExportDir  'ExportDir' names a folder that does
%                                     not exist and could not be created
%
%   EXAMPLE
%     R = demo_gravity();                 % ~10 s, prints the three sections
%     R.report();                         % the toolbox's own verdict
%
%   COST
%     Five days, 16 arcs, three streams, 77 common unknowns and 96 local ones:
%     a few seconds. The SH partials are recomputed per arc rather than
%     cached, which is the honest streaming behaviour and still cheap here.
%
%   See also DEMO_GRAVITY_DATALOADER, DEMO_GRAVITY_DESIGN, DEMO_GRACE,
%   LSQVCE.ADJUSTVCE
%
%   Created 2026-08-21 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

arguments
    opt.MaxOuter (1,1) double {mustBePositive, mustBeInteger} = 4
    opt.Verbose (1,1) logical = false
    opt.Plots (1,1) logical = false
    opt.ExportDir (1,:) char = ''
end

    if ~isempty(opt.ExportDir) && exist(opt.ExportDir, 'dir') ~= 7
        [ok, msg] = mkdir(opt.ExportDir);
        if ~ok
            error('lsqvce:demoGravity:noExportDir', ...
                  'Cannot create ExportDir ''%s'': %s', opt.ExportDir, msg);
        end
    end

    [cfg, L] = demo_gravity_dataLoader();

    % ---- 1. the problem -------------------------------------------------
    section('1. The problem: GRACE-like, physically structured, ragged');
    fprintf('  degree/order   : nmax = %d  ->  %d common unknowns\n', ...
            cfg.nmax, cfg.nc);
    fprintf('  days / arcs    : %d days, %d arcs total (%s per day)\n', ...
            numel(cfg.days), cfg.nArcTot, ...
            strjoin(arrayfun(@(d) sprintf('%d', numel(cfg.arcsOfDay{d})), ...
                             1:numel(cfg.days), 'UniformOutput', false), '/'));
    fprintf('  arc lengths    : %d to %d epochs (every arc its own)\n', ...
            min(cfg.mOfArc(cfg.mOfArc > 0)), max(cfg.mOfArc(:)));
    for s = 1:numel(cfg.sys)
        fprintf('  system %-8s: obsType %-5s  %d local(s)/arc  ->  nl = %d\n', ...
                cfg.sys(s).name, cfg.sys(s).obsType, ...
                cfg.sys(s).nlPerArc, cfg.sys(s).nl);
    end
    fprintf(['  missing arcs   : day 3 has arcs [1 3], day 5 has [2 3 4] --\n' ...
             '                   absent BY ID, the way a real schedule is\n']);

    % ---- run ------------------------------------------------------------
    P = lsqvce.adjustVCE(L, cfg.days, 'MaxOuter', opt.MaxOuter, ...
                         'Verbose', opt.Verbose);
    R = P.groupComponentsBy('obsType').run();

    % ---- 2. the weights --------------------------------------------------
    section('2. The weights VCE had to find');
    fprintf('  %-6s  %10s  %10s  %8s\n', 'stream', 'planted', 'recovered', 'ratio');
    got = zeros(1, numel(R.noise.stream));
    want = zeros(1, numel(R.noise.stream));
    for k = 1:numel(R.noise.stream)
        ot = streamType(R.noise.stream{k});
        got(k)  = R.noise.model{k}.sigma;
        want(k) = cfg.trueSigma.(ot);
        fprintf('  %-6s  %10.4f  %10.4f  %8.3f\n', ot, want(k), got(k), ...
                got(k) / want(k));
    end
    fprintf(['  The ratios matter, not the absolute values: a common factor\n' ...
             '  is absorbed by sigma0 = %.4f.\n'], R.sigma0);

    % ---- 3. the field ----------------------------------------------------
    section('3. The field: recovery read in units of the formal error');
    sig = sqrt(diag(R.sol.Cxx)) * R.sigma0;
    z   = (R.x - cfg.xTrue) ./ sig;
    fprintf('  max |z| over all %d coefficients : %.2f  (5 would be alarming)\n', ...
            cfg.nc, max(abs(z)));
    fprintf('  mean z^2                         : %.3f  (1.00 is ideal)\n', ...
            mean(z.^2));
    kA = cfg.meta.n == 5 & cfg.meta.m == 3 & cfg.meta.cs == 'C';
    fprintf(['  the planted anomaly (5,3,C)      : truth %+.4f, ' ...
             'estimate %+.4f +- %.4f\n'], ...
            cfg.xTrue(kA), R.x(kA), sig(kA));

    [deg, ratio] = degreeVariances(cfg, R.x, sig);
    fprintf('  error degree variance / formal, by degree:\n    ');
    fprintf('n=%d %.2f   ', [deg; ratio]);
    fprintf('\n');

    % ---- plots ------------------------------------------------------------
    if opt.Plots
        f1 = plotDegreeVariances(cfg, deg, ratio);
        f2 = plotCoefficientErrors(cfg, z);
        exportIf(opt.ExportDir, f1, 'demo_gravity_degree_variances');
        exportIf(opt.ExportDir, f2, 'demo_gravity_coefficient_errors');
    end
end

%==========================================================================
function section(title)
%SECTION  Print a section banner. INPUT title 1xK char. OUTPUT none.
    fprintf('\n%s\n%s\n', title, repmat('-', 1, numel(title)));
end

function ot = streamType(streamKey)
%STREAMTYPE  The obsType part of a 'sysID|obsType' stream key.
%   INPUT streamKey 1xK char/string. OUTPUT ot 1xN char.
    ot = char(extractAfter(string(streamKey), '|'));
    if isempty(ot), ot = char(streamKey); end
end

function [deg, ratio] = degreeVariances(cfg, x, sig)
%DEGREEVARIANCES  Realised error power per degree over its formal expectation.
%   INPUT   cfg 1x1 struct; x ncx1 double estimate; sig ncx1 double formal
%   OUTPUT  deg 1xN double the degrees; ratio 1xN double the ratios
    deg   = 2:cfg.nmax;
    ratio = zeros(1, numel(deg));
    err   = x - cfg.xTrue;
    for i = 1:numel(deg)
        k = cfg.meta.n == deg(i);
        ratio(i) = sum(err(k).^2) / sum(sig(k).^2);
    end
end

function f = plotDegreeVariances(cfg, deg, ratio)
%PLOTDEGREEVARIANCES  Realised-over-formal error power per degree.
%   INPUT cfg; deg 1xN; ratio 1xN. OUTPUT f 1x1 matlab.ui.Figure.
    f = figure('Name', 'demo_gravity: error degree variances');
    bar(deg, ratio);
    yline(1, 'r--', 'formal expectation');
    xlabel('spherical-harmonic degree n');
    ylabel('realised / formal error power');
    title(sprintf('Error degree variances (nmax = %d)', cfg.nmax));
    grid on
end

function f = plotCoefficientErrors(cfg, z)
%PLOTCOEFFICIENTERRORS  Normalised error of every coefficient, by degree.
%   INPUT cfg; z ncx1 double. OUTPUT f 1x1 matlab.ui.Figure.
    f = figure('Name', 'demo_gravity: normalised coefficient errors');
    scatter(cfg.meta.n, z, 18, 'filled');
    yline([-3 3], 'r--');
    xlabel('spherical-harmonic degree n');
    ylabel('(estimate - truth) / formal error');
    title('Coefficient recovery in units of the formal error');
    grid on
end

function exportIf(dirName, f, base)
%EXPORTIF  Save a figure as PNG when an export folder was requested.
%   INPUT dirName 1xK char; f figure; base 1xN char. OUTPUT none.
    if isempty(dirName), return; end
    file = fullfile(dirName, [base '.png']);
    exportgraphics(f, file, 'Resolution', 150);
    fprintf('  wrote %s\n', file);
end
