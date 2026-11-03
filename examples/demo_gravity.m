function [R, score] = demo_gravity(opt)
%DEMO_GRAVITY  GRACE-like gravity recovery on a PHYSICALLY STRUCTURED design.
%
%   R = DEMO_GRAVITY()
%   [R, score] = DEMO_GRAVITY()
%   R = DEMO_GRAVITY('MaxPasses', 6, 'Verbose', true)
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
%     2. THE WEIGHTS. The noise declared to the toolbox has the right AR
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
%     'MaxPasses'  1x1 double   pass BUDGET of the estimator; the default
%                               run settles in 4             (default 10)
%     'Verbose'    1x1 logical  print the toolbox's own report (default false)
%     'Plots'      1x1 logical  draw the two summary figures   (default false)
%     'ExportDir'  1xK char     if non-empty, save the figures there as PNG
%                                                              (default '')
%
%   OUTPUT
%     R      1x1 lsqvce.Results  the full result: R.x, R.sigma0, R.sol.Cxx,
%                              R.noise (the recovered stream models), R.meta
%                              (run provenance), R.report() for the verdict.
%     score  1x1 struct        the same self-assessment the demo PRINTS,
%                              returned so it can be checked rather than read
%                              (WP-53). Fields:
%                              .sigma0        1x1 double
%                              .stream        1xS cell    obsType per stream
%                              .sigmaPlanted  1xS double  the planted sigmas
%                              .sigmaRecovered 1xS double what VCE found
%                              .sigmaRatio    1xS double  recovered/planted
%                              .zMax          1x1 double  max |x-xTrue|/sigma
%                              .zMeanSq       1x1 double  mean of z^2
%                              .anomaly       1x1 struct  .true .est .sigma .z
%                              .degree        1xN double  degrees 2..nmax
%                              .degRatio      1xN double  realised/formal
%                                             error power per degree
%
%   ACCURACY
%     MEASURED on the default configuration -- the run this file performs
%     (Estimate="noise", Pool="stream", ArcFactors=false, Noise.Smooth=1,
%     MaxPasses=10) -- on 2026-10-16, MATLAB R2026a, v6.25.2:
%       passes                     4, R.stopReason "converged"
%       sigma0                     1.0000
%       stream sigmas    planted   1.0 / 2.0 / 0.5
%                        recovered 1.1719 / 2.0892 / 0.5080
%                        ratios    1.172 / 1.045 / 1.016
%       max |z| over 77 coefficients   2.32   (5 would be alarming)
%       mean z^2                       0.850  (1.00 ideal; sd = sqrt(2/77)
%                                              = 0.161, so this is 0.9 sd low)
%       planted anomaly (5,3,C)   truth +0.1764, estimate +0.2742 +- 0.0887
%                                 -> 1.10 sigma, comfortably recovered
%       error degree variance / formal, by degree
%         n=2 0.18  n=3 0.98  n=4 0.62  n=5 0.64  n=6 0.87  n=7 1.00  n=8 0.80
%     WHY Noise.Smooth = 1 (WP-145, measured). With the plain 60-bin law the
%     loop did not settle in the four passes this demo used to spend: the
%     zero-frequency bin is what the per-arc locals absorb (redundancy 5-6
%     against a median of 48 per bin), and the update crept toward it at
%     0.977 per pass -- 72 passes to meet the tolerance. Its estimate was
%     statistically as good long before (x within 0.29 sigma_x of the 72-pass
%     fixed point after 4 passes), so the stop was the issue, not the
%     solution. Smoothing each bin with its neighbours costs nothing on these
%     smooth AR spectra (the numbers above match the 60-bin fixed point's:
%     mean z^2 0.920, max |z| 2.28, the same ratios to 1 %) and settles.
%     REPEATABILITY: the noise realisation is seeded, but the Hutchinson
%     probes are NOT (this demo sets no Options.Seed), so the last decimals
%     move -- a second run gave 1.1716 / 2.0892 / 0.5080, max |z| 2.3175,
%     mean z^2 0.8507, sigma0 0.9999, again 4 passes.
%     Reading: the RELATIVE weights come back under REAL column correlation.
%     potD lands 4.5 % high and kbr 1.6 %; potC, the strongly coloured stream
%     (ar [1 -0.6 0.2]), comes back 17 % high -- the bias the shape/level
%     split carries on a short record, the same magnitude the suite bounds at
%     25 % (dev/tests/test_vcePolicy, testRecoversPlantedRelativeWeights).
%     mean z^2 about one sd BELOW 1 means the realised errors are SMALLER than the
%     formal ones: the formal errors are mildly conservative, not
%     over-optimistic, and within the scatter the value is consistent with 1.
%     The n=2 ratio, 0.18 over a degree with only five coefficients, is the
%     least stable of the seven and is quoted as measured, not explained away.
%     These are the numbers of ONE run; they are quoted, not asserted -- the
%     tests assert the structural claims and leave the statistics to the
%     formal-error reading.
%     DECLARATION: measured on the v6 data struct. WP-161 declares the same
%     rows as a lsqvce.TaskDefinition (DEMO_GRAVITY_DATALOADER keeps every
%     seed, draw and the dense local design); not re-measured in the
%     migrating session (no MATLAB there).
%
%   ERRORS
%     lsqvce:demoGravity:noExportDir  'ExportDir' names a folder that does
%                                     not exist and could not be created
%     A declaration mistake is raised by the add-call that makes it
%     (lsqvce:TaskDefinition:*), a mistake inside a Build while the day is
%     read (lsqvce:DayContext:*); the engines' refusals keep their own ids.
%
%   EXAMPLE
%     R = demo_gravity();                 % ~1.6 s, prints the three sections
%     [R, score] = demo_gravity();        % ... and hand the numbers back
%     R.report();                         % the toolbox's own verdict
%
%   COST
%     Five days, 16 arcs, three streams, 77 common unknowns and 96 local ones:
%     about 1.6 s (measured 2026-10-16, R2026a, the four passes the default
%     run settles in). The
%     SH partials are recomputed per arc rather than
%     cached, which is the honest streaming behaviour and still cheap here.
%
%   See also DEMO_GRAVITY_DATALOADER, DEMO_GRAVITY_DESIGN, DEMO_GRACE,
%   LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.OPTIONS
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: the run fits the loader's lsqvce.TaskDefinition. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-150 (F7): the formal error is
%   sqrt(diag(R.sol.Cxx)) -- it was multiplied by sigma0 a second time. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: the run converges --
%   Noise.Smooth 1 and the toolbox's pass budget 10 (it settles in 4; the
%   plain 60-bin law needed 72), the ACCURACY block re-measured. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the ACCURACY block
%   re-measured on the run this file performs today (v5.3, MATLAB R2026a),
%   the documented option name is MaxPasses, and See also drops the deleted
%   engine.adjustVCE. With assistance from Claude (Anthropic), model
%   Claude Opus 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- the run is lsqvce.fit over the v4 data
%   struct (Estimate="noise", Pool="stream", ArcFactors=false); the loaders and the engine call are gone from
%   the example. With assistance from Claude (Anthropic), model Claude
%   Fable 5.
%   Updated 2026-08-21 (CEST, UTC+02:00), WP-53: the printed self-assessment
%   is also RETURNED (second output), and the measured numbers of the default
%   run are pinned in the ACCURACY block above -- the quality debt of the
%   demo's first round, settled by a real run. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Created 2026-08-21 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

arguments
    opt.MaxPasses (1,1) double {mustBePositive, mustBeInteger} = 10
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

    [cfg, task] = demo_gravity_dataLoader();   % the declaration (this folder)

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
    % Noise.Smooth 1 (WP-145, measured): the zero-frequency bin of the 60-bin
    % law is what the per-arc locals absorb (redundancy 5-6 against a median
    % of 48 per bin), so the unsmoothed loop crept toward it at 0.977 per pass
    % and needed 72 passes. Averaging each bin with its two neighbours costs
    % nothing on these smooth AR spectra and settles in 4 passes.
    R = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=opt.MaxPasses, ...
            Noise=struct('Smooth', 1), Verbose=opt.Verbose));

    % ---- 2. the weights --------------------------------------------------
    section('2. The weights VCE had to find');
    fprintf('  %-6s  %10s  %10s  %8s\n', 'stream', 'planted', 'recovered', 'ratio');
    got = zeros(1, numel(R.noise.stream));
    want = zeros(1, numel(R.noise.stream));
    names = cell(1, numel(R.noise.stream));
    for k = 1:numel(R.noise.stream)
        ot = streamType(R.noise.stream{k});
        names{k} = ot;
        got(k)  = R.noise.model{k}.sigma;
        want(k) = cfg.trueSigma.(ot);
        fprintf('  %-6s  %10.4f  %10.4f  %8.3f\n', ot, want(k), got(k), ...
                got(k) / want(k));
    end
    fprintf(['  The ratios matter, not the absolute values: a common factor\n' ...
             '  is absorbed by sigma0 = %.4f.\n'], R.sigma0);

    % ---- 3. the field ----------------------------------------------------
    section('3. The field: recovery read in units of the formal error');
    sig = sqrt(diag(R.sol.Cxx));            % sol.Cxx is sigma0^2-scaled already (WP-150: was x sigma0 again)
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

    % ---- the same numbers, returned rather than only printed (WP-53) -----
    % kA (the anomaly column) was already located for the printed line above.
    score = struct( ...
        'sigma0',          R.sigma0, ...
        'stream',          {names}, ...
        'sigmaPlanted',    want, ...
        'sigmaRecovered',  got, ...
        'sigmaRatio',      got ./ want, ...
        'zMax',            max(abs(z)), ...
        'zMeanSq',         mean(z.^2), ...
        'anomaly',         struct('true', cfg.xTrue(kA), 'est', R.x(kA), ...
                                  'sigma', sig(kA), 'z', z(kA)), ...
        'degree',          deg, ...
        'degRatio',        ratio);

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
