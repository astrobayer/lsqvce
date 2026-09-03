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
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M5, WP-87.
    tests = functiontests(localfunctions);
end

function testKnobCountIsDeliberate(tc)
% categories: contract
    s = lsqvce.Options().toStruct();
    tc.verifyEqual(numel(fieldnames(s)), 24, ...
        '13 knobs + 5 groups + 6 interim (until M5 5b) -- update on purpose');
    tc.verifyEqual(sort(fieldnames(s)).', sort({'Estimate','Pool','ArcFactors','Robust', ...
        'MaxPasses','MinPasses','Tolerance','Backend','Probes','Seed','CacheDays','Verbose','LogFile', ...
        'Noise','RobustOpts','Diagnostics','Parallel','GaussNewton', ...
        'Reg','RegLambda','RegMatrix','Sparse','DatumB','Sigma0'}));
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

function testGroupMergeOntoDefaults(tc)
% categories: contract
    o = lsqvce.Options(Noise=struct('Grid', 2160));
    tc.verifyEqual(o.Noise.Grid, 2160);
    tc.verifyEqual(o.Noise.Whitener, "ar", 'untouched fields keep their defaults');
    tc.verifyEqual(o.Noise.SkipThreshold, 0.2);
    o.Parallel.Mode = "parfor";                       % field assignment
    tc.verifyEqual(o.Parallel.Mode, "parfor");
    tc.verifyEqual(o.Parallel.Threshold, 8);
    tc.verifyError(@() lsqvce.Options(Noise=struct('Lags', 5)), 'lsqvce:Options:badGroupField');
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
