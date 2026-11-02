function R = demo_gravfire(bundleFile, opt)
%DEMO_GRAVFIRE  A GRAVFIRE bundle (commonData / myOptions / systemsData) through fit.
%   R = demo_gravfire(bundleFile)
%   R = demo_gravfire(bundleFile, Estimate="noise", Passes=3, Verbose=1)
%
%   The reading example the R5 track was written for: a GRACE two-day bundle
%   (Stokes coefficients as the common parameters, per-day ICs / accelerometer
%   parameters / empiricals as the locals, kinematic positions of both
%   satellites and the K-band range-rate as three streams) goes through
%   template/gravfire_data.m into a lsqvce.TaskDefinition and through
%   lsqvce.fit. The design function is NOT part of the toolbox (the
%   DesignFcn of the guide is yours); the template's 'DryRun' stands in for
%   it with random design matrices, so this demo exercises the PIPELINE --
%   the declaration, the day-scoped locals, the per-epoch covariance of the
%   3-axis positions, the joint estimator with Options.OnPass -- not gravity
%   recovery.
%
%   The bundle is your private file (testdata.mat, ~84 MB) and is not shipped:
%   without it this demo prints what it would do and returns [] -- it never
%   fails a smoke run.
%
%   INPUT   bundleFile  1x1 string  path to testdata.mat ("" = look next to the toolbox)
%   NAME-VALUE
%     Estimate  1x1 string  "none" | "noise"   (default "none": the pipeline check -- one solve)
%     Passes    1x1 double  MaxPasses for Estimate "noise" (default 3)
%     Verbose   1x1 double  0..2 (default 1)
%   OUTPUT  R  1x1 lsqvce.Results | []   [] when there is no bundle at bundleFile
%
%   ACCURACY
%     Not applicable: a pipeline check on random design matrices; its numbers
%     mean nothing (the template's dry run).
%
%   ERRORS  none of its own; the template's gravfire:* and the declaration's
%           lsqvce:TaskDefinition:* pass through (help gravfire_data).
%
%   See also TEMPLATE/GRAVFIRE_DATA, LSQVCE.TASKDEFINITION, LSQVCE.FIT,
%            dev/GRAVFIRE_INTERFACE.md
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: the bundle becomes a
%   lsqvce.TaskDefinition (gravfire_data's task) and fit reads the task. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the v5 data struct, the
%   machine-specific timing dropped from the Estimate line, and Verbose 0..2
%   (level 3 was merged into 2). With assistance from Claude (Anthropic),
%   model Claude Opus 5.
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M6, WP-88 (gracefoReading).
%   Updated 2026-10-29 (CET, UTC+01:00), WP-159: the declaration is the data
%   struct of the API reference (it was called "v5"; since v6 a day is one
%   cell per arc). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
arguments
    bundleFile (1,1) string = ""
    opt.Estimate (1,1) string {mustBeMember(opt.Estimate, ["none","noise"])} = "none"
    opt.Passes (1,1) double {mustBeInteger, mustBePositive} = 3
    opt.Verbose (1,1) double = 1
end
    root = fileparts(fileparts(mfilename('fullpath')));
    if strlength(bundleFile) == 0
        bundleFile = fullfile(root, '..', 'testdata.mat');
    end
    if ~isfile(bundleFile)
        fprintf(['demo_gravfire: no bundle at %s.\n' ...
                 '  Point it at your GRAVFIRE testdata.mat (commonData, myOptions, systemsData).\n' ...
                 '  It would: copy template/gravfire_data.m into your project, declare the\n' ...
                 '  task (DryRun = random design), and run lsqvce.fit with Estimate="%s".\n'], ...
                 bundleFile, opt.Estimate);
        R = [];
        return
    end
    addpath(fullfile(root, 'template'));               % gravfire_data.m (copy it into your project)
    S = load(bundleFile);
    task = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
    watch = @(ctx) fprintf('  pass %d: max|alpha2-1| = %.3g   levels = %s\n', ...
                           ctx.pass, ctx.maxDev, mat2str(ctx.sigma, 3));
    if opt.Estimate == "none"
        R = lsqvce.fit(task, Verbose=opt.Verbose);
    else
        w = warning('off', 'lsqvce:adjustVCEjoint:notConverged');  c = onCleanup(@() warning(w)); %#ok<NASGU>
        R = lsqvce.fit(task, Estimate="noise", Pool="obsType", ...
                       Noise=struct('Grid', 2160), MaxPasses=opt.Passes, CacheDays=Inf, ...
                       Verbose=opt.Verbose, OnPass=watch);
    end
    if opt.Verbose > 0, R.report; end
end
