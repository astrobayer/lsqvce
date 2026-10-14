function tests = test_options
%TEST_OPTIONS  lsqvce.Options v5: 13 knobs + 5 groups, typed, one mechanism.
%   tests = test_options
%
%   R5 M5 (WP-87, decision R5-D2): the option space was cut from 62 flat knobs
%   to 13 top-level properties and five typed groups; Name=Value on
%   lsqvce.fit is the same mechanism as the object. Pinned here: the count
%   (deliberately), typed refusals with MATLAB ids at assignment, the group
%   merge (a partial struct onto the defaults; unknown fields refused with
%   the lsqvce id), value semantics, the engine argument lists, and the
%   display.
%
%   See also LSQVCE.OPTIONS, LSQVCE.FIT, TEST_FIT
%
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135: Noise.Whitener validated in
%   any case. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: fit's pairs override opts
%   field by field. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M5, WP-87.
    tests = functiontests(localfunctions);
end

function testFitPairsOverrideOptsFieldByField(tc)
% categories: contract
    % WP-131: fit's help says the Name=Value pairs override OPTS. A partial
    % group struct went through the group setter, which merges onto the
    % DEFAULTS: fit(data, Options(Noise=struct('Whitener',"ar")),
    % Noise=struct('Grid',16)) ran with the default whitener. fit now merges
    % the pair onto opts' group first; the Options object itself keeps its
    % documented semantics (a partial struct onto the defaults).
    rs = RandStream('mt19937ar', 'Seed', 131);
    A = randn(rs, 30, 3);  l = A*[1; 2; 3] + 0.1*randn(rs, 30, 1);
    o = lsqvce.Options(Noise=struct('Whitener', "ar"), Verbose=false);
    R = lsqvce.fit(struct('A', A, 'l', l), o, Noise=struct('Grid', 16));
    tc.verifyEqual(string(R.meta.options.Noise.Whitener), "ar", 'the field opts set survives');
    tc.verifyEqual(R.meta.options.Noise.Grid, 16, 'the field the pair set wins');
    o2 = lsqvce.Options(Noise=struct('Whitener', "ar"));
    o2.Noise = struct('Grid', 16);
    tc.verifyEqual(string(o2.Noise.Whitener), "toeplitz", 'the object keeps merging onto the defaults');
end

function testKnobCountIsDeliberate(tc)
% categories: contract
    s = lsqvce.Options().toStruct();
    tc.verifyEqual(numel(fieldnames(s)), 17, ...
        ['12 knobs + 5 groups -- update on purpose. Reg/RegLambda/RegMatrix went ' ...
         'in R5-D8 round 5 (regularisation is DECLARED); MinPasses and LogFile ' ...
         'went in R6 WP-101, where they were found to have no consumer anywhere ' ...
         'in the package -- help line, property, disp table, and nothing else']);
    tc.verifyEqual(sort(fieldnames(s)).', sort({'Estimate','Pool','ArcFactors','Robust', ...
        'MaxPasses','Tolerance','Backend','Probes','Seed','CacheDays','Verbose','OnPass', ...
        'Noise','RobustOpts','Diagnostics','Parallel','GaussNewton'}));
end

function testTypedRefusalsAtAssignment(tc)
% categories: contract, robustness
    tc.verifyError(@() lsqvce.Options(Estimate="vce"),   'MATLAB:validators:mustBeMember');
    tc.verifyError(@() lsqvce.Options(Robust="irls"),    'MATLAB:validators:mustBeMember');
    tc.verifyError(@() lsqvce.Options(MaxPasses=2.5),    'MATLAB:validators:mustBeInteger');
    tc.verifyError(@() lsqvce.Options(Tolerance=0),      'MATLAB:validators:mustBePositive');
    tc.verifyError(@() lsqvce.Options(Backend="lu"),     'MATLAB:validators:mustBeMember');
    tc.verifyError(@() lsqvce.Options(Verbose=7),        'lsqvce:Options:badVerbose');
    tc.verifyError(@() lsqvce.Options(Noise=struct('Grid', "lots")), 'lsqvce:Options:badGrid');
end

function testWhitenerIsValidatedInAnyCase(tc)
% categories: contract, robustness
    % WP-135 (review 3, lane 4): Noise.Whitener was forwarded raw -- "AR"
    % worked under Estimate "none" and died mid-run under "noise"
    % (jointModelOpt's mustBeMember), "foo" ran until an engine refused it.
    tc.verifyError(@() lsqvce.Options(Noise=struct('Whitener', "foo")), 'lsqvce:Options:badWhitener');
    tc.verifyError(@() lsqvce.Options(Noise=struct('Whitener', 3)), 'lsqvce:Options:badWhitener');
    o = lsqvce.Options(Noise=struct('Whitener', "AR"));
    a = o.engineArgs('adjustBlock');
    tc.verifyEqual(a.nv{find(strcmp(a.nv(1:2:end), 'PsdMode')) * 2}, "ar");
end

function testGroupMergeOntoDefaults(tc)
% categories: contract
    o = lsqvce.Options(Noise=struct('Grid', 2160));
    tc.verifyEqual(o.Noise.Grid, 2160);
    tc.verifyEqual(o.Noise.Whitener, "toeplitz", 'untouched fields keep their defaults');
    tc.verifyEqual(o.Noise.SkipThreshold, 0.2);
    o.Parallel.Mode = "parfor";                       % field assignment
    tc.verifyEqual(o.Parallel.Mode, "parfor");
    tc.verifyEqual(o.Parallel.Threshold, 8);
    tc.verifyError(@() lsqvce.Options(Noise=struct('Lags', 5)), 'lsqvce:Options:badGroupField');
    % WP-118 part 2 (user decision 2026-10-10): Noise.MeanType is retired -- the
    % arc factors are absolute since 6.1.6, gamma-bar no longer scales them
    tc.verifyError(@() lsqvce.Options(Noise=struct('MeanType', "plain")), 'lsqvce:Options:badGroupField', ...
        'the retired MeanType is refused by name, not silently ignored');
    tc.verifyFalse(isfield(lsqvce.Options().Noise, 'MeanType'));
    tc.verifyError(@() lsqvce.Options(Diagnostics=42), 'MATLAB:validation:UnableToConvert', 'a non-struct is refused by the type');
end

function testValueSemanticsAndNameOrderInvariance(tc)
% categories: robustness
    o1 = lsqvce.Options(Estimate="noise", Backend="normal", MaxPasses=3);
    a1 = o1.engineArgs("adjustVCEjoint");
    o1.MaxPasses = 9;
    % tol: exact -- value semantics, no floating point involved
    tc.verifyEqual(a1.nv{find(strcmp(a1.nv, 'MaxIter')) + 1}, 3);
    o2 = lsqvce.Options(MaxPasses=3, Backend="normal", Estimate="noise");
    o3 = lsqvce.Options(Estimate="noise", Backend="normal", MaxPasses=3);
    tc.verifyTrue(isequal(o2, o3), 'Name-Value order must not matter');
    tc.verifyEqual(o2.engineArgs("adjustBlock"), o3.engineArgs("adjustBlock"));
end

function testEngineArgumentLists(tc)
% categories: contract, reference
    o = lsqvce.Options(Estimate="noise", Pool="obsType", MaxPasses=4, Tolerance=1e-3, Probes=50, Seed=7);
    a = o.engineArgs("adjustVCEjoint");
    nv = a.nv;
    get = @(name) nv{find(strcmp(nv, name)) + 1};
    tc.verifyEqual(get('NumLags'), 60);
    tc.verifyEqual(get('MaxIter'), 4);
    tc.verifyEqual(get('Tol'), 1e-3);
    tc.verifyEqual(get('Probes'), 50);
    tc.verifyEqual(get('Seed'), 7);
    tc.verifyEqual(get('Pool'), "obsType");
    tc.verifyEqual(get('SkipThreshold'), 0.2);
    r = lsqvce.Options(Robust="igg3", RobustOpts=struct('MaxIter', 12)).engineArgs("adjustRobust");
    tc.verifyEqual(r.nv{find(strcmp(r.nv, 'Influence')) + 1}, "igg3");
    tc.verifyEqual(r.nv{find(strcmp(r.nv, 'MaxIter')) + 1}, 12);
    p = lsqvce.Options(Parallel=struct('Mode', "parfor", 'Threshold', 2)).engineArgs("adjustBlock");
    tc.verifyEqual(p.nv{find(strcmp(p.nv, 'Parallel')) + 1}, "parfor");
    tc.verifyEqual(p.nv{find(strcmp(p.nv, 'ParallelThreshold')) + 1}, 2);
    % WP-125: Noise.Whitener reaches the declared psd laws of a single solve
    % and of the stand-alone IRLS, not only the joint estimator's law
    tc.verifyEqual(p.nv{find(strcmp(p.nv, 'PsdMode')) + 1}, "toeplitz");
    w = lsqvce.Options(Robust="huber", Noise=struct('Whitener', "ar")).engineArgs("adjustRobust");
    tc.verifyEqual(w.nv{find(strcmp(w.nv, 'PsdMode')) + 1}, "ar");
    id = '';
    try, o.engineArgs("nope"); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badTarget');
end

function testNameValueOnFitIsTheSameMechanism(tc)
% categories: contract, reference
    % fit(data, Name=Value) == fit(data, Options(Name=Value)): the resolved
    % options in the provenance are identical, and so is the estimate
    rs = RandStream('mt19937ar', 'Seed', 5);
    A = randn(rs, 40, 3);  l = A*[1; -2; 0.5] + 0.01*randn(rs, 40, 1);
    data = struct('A', A, 'l', l);
    R1 = lsqvce.fit(data, lsqvce.Options(Backend="normal", Verbose=false));
    R2 = lsqvce.fit(data, Backend="normal", Verbose=false);
    R3 = lsqvce.fit(data, lsqvce.Options(Verbose=false), Backend="normal");
    % tol: exact -- one code path, one option object after resolution
    tc.verifyEqual(R2.x, R1.x);
    tc.verifyEqual(R3.x, R1.x);
    tc.verifyEqual(R2.meta.options.Backend, "normal");
    tc.verifyEqual(rmfield(R2.meta.options, 'kind'), rmfield(R1.meta.options, 'kind'));
end

function testDisplayMarksDefaults(tc)
% categories: contract | tags: diagnostics
    o = lsqvce.Options();
    o.MaxPasses = 4;
    txt = evalc('disp(o)');
    tc.verifyTrue(contains(txt, 'lsqvce.Options'));
    tc.verifyTrue(contains(txt, '     MaxPasses'), 'a changed value is listed');
    tc.verifyFalse(contains(txt, '   * MaxPasses'), 'a changed value must NOT carry the star');
    tc.verifyTrue(contains(txt, '   * Backend'), 'a default carries the star');
end
