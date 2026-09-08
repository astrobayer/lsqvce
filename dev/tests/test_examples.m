function tests = test_examples
%TEST_EXAMPLES  Every script in examples/ runs, once, in isolation.
%   tests = test_examples
%
%   WHY THIS SUITE EXISTS. The postflight cited in every pull request of review
%   R6 ends "examples N/N", and until this round that N came from a harness
%   that lived in a scratch directory and was rebuilt by hand each round: the
%   claim had no authority in the tree, and 8 of the 15 files in examples/ were
%   executed by NOTHING here -- bench_qr_vs_normal, demo_combine, demo_grace,
%   demo_gravfire, demo_lsqvce, demo_messaging, demo_sensors, demo_workflow.
%   A demo nobody runs rots, and R6 spent two rounds proving it: WP-95 found 54
%   confirmed defects in what the examples SAID, and every one of them had
%   survived a green gate.
%
%   WHAT IT ASSERTS, and the two halves are not the same claim:
%     * every example RUNS -- one call each, in its own workspace, figures off;
%     * every example is IN THE TABLE -- a call table nobody updates is worse
%       than no table, so adding examples/demo_new.m without a row FAILS here,
%       naming the file.
%
%   WHAT IT DOES NOT ASSERT. Nothing about the numbers a demo prints. Those
%   belong to the demo's own guarding test (demo_fit in test_fit, demo_gravity
%   in test_physOracle, demo_priors in test_fit, demo_noisePolicy in
%   test_vcePolicy), where the planted truth is. This suite is a smoke test and
%   says so.
%
%   ISOLATION MATTERS, measured: the first version of this harness ran every
%   example in ONE shared workspace and clobbered itself -- a script's own
%   variables overwrote the loop counter and the run died with "Unrecognized
%   function or variable 'nm'". Each call therefore runs inside its own
%   function (runOne / runScript) and figures are closed between them.
%
%   See also RUN_ALL_LSQVCE, CI_GATE, TEST_PHYSORACLE
%
%   Created 2026-09-08 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5 -- R6 WP-99: the example harness had lived
%   in a scratch directory all review long (HANDOVER 12.23).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    tc.TestData.root  = root;
    tc.TestData.exDir = fullfile(root, 'examples');
    addpath(root);  addpath(tc.TestData.exDir);  addpath(fullfile(root, 'template'));
    vis = get(groot, 'DefaultFigureVisible');
    set(groot, 'DefaultFigureVisible', 'off');
    tc.addTeardown(@() set(groot, 'DefaultFigureVisible', vis));
    tc.addTeardown(@() close('all', 'force'));
end

function testEveryExampleIsInTheCallTable(tc)
% categories: contract
    % The half that keeps the other half honest. exampleCalls is a hand-written
    % table because the calls are not uniform -- two examples take required
    % arguments and two are scripts -- so the table can only be trusted if
    % nothing can be added to examples/ without appearing in it.
    onDisk = dir(fullfile(tc.TestData.exDir, '*.m'));
    onDisk = sort({onDisk.name});
    onDisk = cellfun(@(f) f(1:end-2), onDisk, 'UniformOutput', false);
    inTable = sort(exampleCalls(tc.TestData.exDir, 'names'));
    tc.verifyEqual(inTable, onDisk, ...
        sprintf(['examples/ and the call table disagree. Missing from the table: %s. ' ...
                 'In the table but not on disk: %s.'], ...
                 strjoin([setdiff(onDisk, inTable), {'(none)'}], ', '), ...
                 strjoin([setdiff(inTable, onDisk), {'(none)'}], ', ')));
end

function testEveryExampleRuns(tc)
% categories: robustness
    % One call each, every failure collected before any is reported: a run that
    % stops at the first broken demo hides the other fourteen.
    [names, calls] = exampleCalls(tc.TestData.exDir, 'both');
    bad = {};
    for k = 1:numel(names)
        try
            runOne(calls{k});
        catch ME
            bad{end+1} = sprintf('%s (%s: %s)', names{k}, ME.identifier, ...
                                 strrep(ME.message, newline, ' ')); %#ok<AGROW>
        end
        close('all', 'force');
    end
    tc.verifyEmpty(bad, sprintf('examples that failed to run: %s', strjoin(bad, ' | ')));
end

%==========================================================================
function varargout = exampleCalls(exDir, what)
%EXAMPLECALLS  THE table: one row per file in examples/, with the call that
%   exercises it. Two examples take required arguments and two are scripts, so
%   this cannot be derived from the folder listing -- which is why the census
%   test above refuses to let it drift. Arguments are chosen for SPEED: this
%   suite is a smoke test inside a gate that already takes minutes, so every
%   knob that shortens a run without changing what is exercised is turned down.
    orbit = struct('u', linspace(0, 4*pi, 400).', 'lam0', linspace(0, 0.3, 400).', ...
                   'incl', 89*pi/180, 'r', 1.07);
    T = { ...
      'bench_qr_vs_normal',      @() bench_qr_vs_normal('u', [10 20], 'n', 2e4, 'nk', 2000); ...
      'demo_combine',            @() demo_combine(); ...
      'demo_fit',                @() demo_fit(); ...
      'demo_grace',              @() demo_grace('verbose', false); ...
      'demo_grace_dataLoader',   @() demo_grace_dataLoader(); ...
      'demo_gravfire',           @() demo_gravfire(""); ...
      'demo_gravity',            @() demo_gravity(); ...
      'demo_gravity_dataLoader', @() demo_gravity_dataLoader(); ...
      'demo_gravity_design',     @() demo_gravity_design(4, orbit); ...
      'demo_lsqvce',             @() runScript('demo_lsqvce'); ...
      'demo_messaging',          @() demo_messaging(); ...
      'demo_noisePolicy',        @() demo_noisePolicy(); ...
      'demo_priors',             @() demo_priors(); ...
      'demo_sensors',            @() demo_sensors(); ...
      'demo_workflow',           @() runScript('demo_workflow')};
    if strcmp(what, 'names')
        varargout = {sort(T(:, 1).')};
    else
        varargout = {T(:, 1).', T(:, 2).'};
    end
    assert(~isempty(exDir), 'the examples folder must be known');
end

function runOne(fh)
%RUNONE  One example, silently, in a workspace of its own.
    evalc('fh()');
end

function runScript(name)
%RUNSCRIPT  A SCRIPT example: it needs a workspace to write its variables into,
%   and this function is it. Sharing one would let a script overwrite the
%   harness's own variables -- measured, it did.
    evalc(name);
end
