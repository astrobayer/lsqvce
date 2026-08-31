function tests = test_options
%TEST_OPTIONS  Tests for lsqvce.Options, the one settings object (R4 WP-55).
%
%     lsqvce.Options   typed knobs with the engines' defaults; engineArgs maps
%                      the object onto one existing driver's inputs
%
%   The claim under test is EXACTNESS of the mapping: a default Options must
%   reproduce every driver's own defaults byte for byte (no behaviour change,
%   HANDOVER §11.6 WP-55), and every bad value must be rejected at the
%   signature (MATLAB ids) or, for cross-property rules, in the body with an
%   lsqvce id (the b1-r3 rule).
%
%   Run: runtests('tests/test_options')
%
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L1, WP-55.
    tests = functiontests(localfunctions);
end

%% ===================================================== reference: defaults
function testDefaultsMatchAdjustVCEEngine(tc)
% categories: reference, contract
    % The adjustVCE constructor is the authority for the four configurator
    % groups (adjustVCE.m, "obj.shapeOpt = struct(...)" and neighbours). A
    % fresh engine and a default Options must agree field for field, in the
    % engine's own field order and value classes -- isequal on the structs.
    o = lsqvce.Options();
    a = o.engineArgs("adjustVCE");
    sys = struct('sysID', 1, 'nl', 1, 'nc', 2, 'loader', @(varargin) struct());
    P = lsqvce.adjustVCE(sys, 1);                 % construction reads no data
    % tol: exact -- these are the engine's literal defaults, not measurements
    tc.verifyEqual(a.shape, P.shapeOpt, 'stochastic-model group must equal the engine default');
    tc.verifyEqual(a.scale, P.scaleOpt, 'variance-component group must equal the engine default');
    tc.verifyEqual(a.outer, P.outerOpt, 'outer-loop group must equal the engine default');
    tc.verifyEqual(a.joint, P.jointOpt, 'joint group must equal the engine default');
    tc.verifyEqual(a.mode, P.vceMode, 'default Model "gls"/"vce" must map to ''separate''');
    tc.verifyEqual(a.axes, struct('pos', 3), 'components-per-epoch default must be pos -> 3');
    tc.verifyEqual(P.axesMap('pos'), 3);
    tc.verifyEqual(struct('Backend', char(o.Backend), 'Reg', char(o.Reg), ...
                          'Parallel', char(o.Parallel)), P.blockOpt, ...
        'solver group must equal the engine''s blockOpt');
end

function testDriverDefaultsArePinned(tc)
% categories: reference
    % Literal defaults READ FROM THE DRIVER SIGNATURES on 2026-08-31:
    % adjustBlock.m arguments block, adjust.m:90-95, robust.m:115-123,
    % solveNonlinear.m:121-128, adjustVCEjoint.m:173-184, diagnose.m:100-111.
    % A driver default that moves breaks this test ON PURPOSE: Options must
    % then be updated in the same edit (one authority per default).
    % ONE pinned deviation: adjustBlock's own Verbose default is FALSE;
    % Options.Verbose defaults to TRUE (the adjust/adjustVCE convention), so
    % the adjustBlock pin below carries true -- documented in the class help.
    o = lsqvce.Options();
    % tol: exact -- literal pins
    tc.verifyEqual(nv2struct(o.engineArgs("adjustBlock").nv), struct( ...
        'Backend', "qr", 'Reg', "none", 'RegLambda', 0, 'RegMatrix', [], ...
        'Parallel', "serial", 'Verbose', true, 'ShowSeconds', 0, 'ExportDir', ""));
    tc.verifyEqual(nv2struct(o.engineArgs("adjustVCEjoint").nv), struct( ...
        'NumLags', 60, 'MaxIter', 10, 'Tol', 1e-2, 'Probes', 100, 'Seed', [], ...
        'Norm', "ortho", 'MeanType', "robust", 'StepSchedule', "none", ...
        'AvgTail', 0, 'Verbose', true));
    d = struct('Alpha', 1e-3, 'Bonferroni', true, 'Beta', 0.20, 'RMin', 1e-6, ...
               'RWeak', 0.10, 'Test', "auto", 'RobustCheck', true, 'Store', "all", ...
               'NWorst', 1000, 'Verify', true, 'Progress', false);
    tc.verifyEqual(nv2struct(o.engineArgs("diagnose").nv), d);
    b = nv2struct(o.engineArgs("adjust").nv);
    tc.verifyEqual(rmfield(b, fieldnames(d)), struct('Method', 'qr', 'Sparse', false, ...
        'DatumB', [], 'Sigma0', 1, 'Verbose', true));
    r = nv2struct(o.engineArgs("robust").nv);
    tc.verifyEqual(rmfield(r, fieldnames(d)), struct('Influence', 'huber', 'Tuning', [], ...
        'MaxIter', 30, 'Tol', 1e-8, 'DatumB', [], 'Method', 'qr', 'Sigma0', 1, 'Verbose', true));
    n = nv2struct(o.engineArgs("solveNonlinear").nv);
    tc.verifyEqual(rmfield(n, fieldnames(d)), struct('MaxIter', 30, 'TolX', 1e-3, ...
        'TolF', 1e-10, 'Lambda0', 1e-6, 'DatumB', [], 'Method', 'qr', 'Sigma0', 1, 'Verbose', true));
    % the diagnose pass-through rides on every batch driver (replaces KeepUnmatched)
    for f = fieldnames(d).'
        tc.verifyEqual(b.(f{1}), d.(f{1}));
        tc.verifyEqual(r.(f{1}), d.(f{1}));
        tc.verifyEqual(n.(f{1}), d.(f{1}));
    end
end

%% ============================================================== contract
function testSignatureRejectsBadValuesAndUnknownNames(tc)
% categories: contract
    % WP-39 lever: enum VALUES are case-sensitive lowercase (mustBeMember),
    % typed constraints fire at the signature, unknown names are rejected.
    % try/catch identifier pattern (never the verifyError output form).
    ok = false;
    try
        lsqvce.Options(Model="GLS");
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'Model value "GLS" must be rejected (case-sensitive mustBeMember)');
    ok = false;
    try
        lsqvce.Options(RegLambda=-1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'negative RegLambda must be rejected by the signature');
    ok = false;
    try
        lsqvce.Options(NoSuchOption=1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'an unknown option name must be rejected');
    ok = false;
    try
        lsqvce.Options(VceMaxIter=2.5);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'a non-integer iteration count must be rejected');
end

function testLaterAssignmentIsValidated(tc)
% categories: contract
    % Property validation must hold AFTER construction too -- that is the
    % reason Options is a class and not a struct.
    o = lsqvce.Options();
    ok = false;
    try
        o.Backend = "QR";
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'assigning an invalid enum value later must be rejected');
    ok = false;
    try
        o.WhiteAlpha = 1.5;
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'WhiteAlpha >= 1 must be rejected');
    o.MaxOuter = 4;                                  % a valid change sticks
    tc.verifyEqual(o.engineArgs("adjustVCE").outer.MaxOuter, 4);
end

function testCrossPropertyRulesUseLsqvceIds(tc)
% categories: contract
    % b1-r3 rule: what the property validators cannot express is checked in
    % the body with lsqvce ids.
    id = '';
    try, lsqvce.Options(GroupBy=42); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badGroupBy');
    id = '';
    try, lsqvce.Options(NumLags="lots"); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badNumLags');
    id = '';
    try, lsqvce.Options(NumLags=1); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badNumLags');
    id = '';
    try, lsqvce.Options(Verbose=7); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badVerbose');
    o = lsqvce.Options();
    id = '';
    try, o.engineArgs("nope"); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badTarget');
    o = lsqvce.Options(NumLags="auto");             % legal for the pooled engine
    tc.verifyEqual(o.engineArgs("adjustVCE").joint.NumLags, 'auto');
    id = '';
    try, o.engineArgs("adjustVCEjoint"); catch e, id = e.identifier; end
    tc.verifyEqual(id, 'lsqvce:Options:badNumLags', ...
        'adjustVCEjoint needs a numeric NumLags; "auto" must be refused with the lsqvce id');
    % a function-handle GroupBy is accepted and passed through untouched
    g = @(sysID, obsType, arc, day) sprintf('%s|%d', obsType, arc);
    o = lsqvce.Options(GroupBy=g);
    a = o.engineArgs("adjustVCE");
    k = find(strcmp(a.nv, 'GroupBy'), 1);
    tc.verifyTrue(isa(a.nv{k+1}, 'function_handle'));
end

function testModelSelectsTheVceMode(tc)
% categories: contract
    % "pooled" is adjustVCE's 'joint' mode; "vce" and "gls" leave it 'separate'.
    tc.verifyEqual(lsqvce.Options(Model="pooled").engineArgs("adjustVCE").mode, 'joint');
    tc.verifyEqual(lsqvce.Options(Model="vce").engineArgs("adjustVCE").mode, 'separate');
    tc.verifyEqual(lsqvce.Options(Model="gls").engineArgs("adjustVCE").mode, 'separate');
    nv = lsqvce.Options(Model="pooled").engineArgs("adjustVCE").nv;
    k = find(strcmp(nv, 'VceMode'), 1);
    tc.verifyEqual(nv{k+1}, "joint");
end

function testDispMarksEngineDefaults(tc)
% categories: contract | tags: diagnostics
    % The display is a beginner's map of the option space: one line per
    % knob, grouped, with a star on every value that IS the engine default.
    o = lsqvce.Options();
    o.MaxOuter = 4;
    txt = evalc('disp(o)');
    tc.verifyTrue(contains(txt, 'lsqvce.Options'));
    tc.verifyTrue(contains(txt, 'Outer loop'));
    tc.verifyTrue(contains(txt, '   * Model'), 'a default value must carry the star');
    tc.verifyTrue(contains(txt, '     MaxOuter'), 'a changed value is listed');
    tc.verifyFalse(contains(txt, '   * MaxOuter'), 'a changed value must NOT carry the star');
end

%% ============================================================ robustness
function testValueSemanticsAndNameOrderInvariance(tc)
% categories: robustness
    % Options is a VALUE class: args returned earlier do not change when the
    % object is edited later, and the order of Name-Value pairs is irrelevant.
    o1 = lsqvce.Options(Model="vce", Backend="normal", MaxOuter=3);
    a1 = o1.engineArgs("adjustVCE");
    o1.MaxOuter = 9;
    % tol: exact -- value semantics, no floating point involved
    tc.verifyEqual(a1.outer.MaxOuter, 3);
    o2 = lsqvce.Options(MaxOuter=3, Backend="normal", Model="vce");
    o3 = lsqvce.Options(Model="vce", Backend="normal", MaxOuter=3);
    tc.verifyTrue(isequal(o2, o3), 'Name-Value order must not matter');
    tc.verifyEqual(o2.engineArgs("adjustBlock"), o3.engineArgs("adjustBlock"));
end

function testToStructCoversEveryPropertyOnce(tc)
% categories: contract | tags: provenance
    % R4 WP-56b: toStruct() is what fit records in Results.meta.options; it
    % must list every public property exactly once and nothing else.
    o = lsqvce.Options(Model="vce", MaxOuter=4);
    s = o.toStruct();
    tc.verifyEqual(s.Model, "vce");
    tc.verifyEqual(s.MaxOuter, 4);
    txt = evalc('disp(o)');
    for f = fieldnames(s).'
        tc.verifyTrue(contains(txt, [' ' f{1} ' ']) || contains(txt, [' ' f{1} sprintf('\n')]) ...
            || ~isempty(regexp(txt, ['\s' f{1} '\s'], 'once')), ...
            sprintf('property %s must appear in the grouped display', f{1}));
    end
    tc.verifyEqual(numel(fieldnames(s)), numel(unique(fieldnames(s))));
    tc.verifyEqual(numel(fieldnames(s)), 57, 'the option space has 57 knobs (update on purpose)');
end

%% ================================================================ helpers
function s = nv2struct(nv)
% Turn a Name-Value cell into a struct (keeps the driver's value classes).
    s = struct();
    for k = 1:2:numel(nv)
        s.(nv{k}) = nv{k+1};
    end
end
