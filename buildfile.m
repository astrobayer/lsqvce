function plan = buildfile
%BUILDFILE  MATLAB buildtool plan for lsqvce: audit, test, package (R4 WP-65).
%
%   buildtool audit      python3 dev/docs_src/audit_tree.py (exit 0 = clean)
%   buildtool test       the full MATLAB suite through dev/tests/run_all_lsqvce
%   buildtool package    lsqvce_<version>.mltbx in ./release, programmatic
%                        (matlab.addons.toolbox.ToolboxOptions, no .prj)
%   buildtool            audit, then test
%
%   Verified by the user's local postflight (buildtool is not part of the CI
%   gate, which runs dev/tests/ci_gate.m directly).
%
%   Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L5, WP-65.
    plan = buildplan(localfunctions);
    plan("test").Dependencies = "audit";
    plan.DefaultTasks = "test";
end

function auditTask(~)
% Run the tree audit; a finding fails the build.
    [status, out] = system('python3 dev/docs_src/audit_tree.py');
    disp(out);
    assert(status == 0, 'lsqvce:build:audit', 'audit_tree.py reported findings.');
end

function testTask(~)
% Run the full suite (the CI gate's runner).
    % L8 (WP-69): temporary path, restored when the task returns.
    pathBefore  = path();
    restorePath = onCleanup(@() path(pathBefore)); %#ok<NASGU>
    addpath(fullfile(pwd, 'dev', 'tests'));
    run_all_lsqvce();
end

function packageTask(~)
% Package the toolbox as an .mltbx into ./release.
    root = pwd;
    txt  = fileread(fullfile(root, '+lsqvce', 'Contents.m'));
    ver  = regexp(txt, '^%\s*Version\s+(\S+)', 'tokens', 'once', 'lineanchors');
    if isempty(ver), ver = {'0.0'}; end
    files = [dir(fullfile(root, '+lsqvce', '**', '*.m')); ...
             dir(fullfile(root, 'html', '*')); ...
             dir(fullfile(root, 'examples', '*.m'))];
    paths = string(arrayfun(@(f) fullfile(f.folder, f.name), files(~[files.isdir]), ...
                            'UniformOutput', false));
    paths = [paths; fullfile(root, "info.xml"); fullfile(root, "helptoc.xml"); ...
             fullfile(root, "README.md"); fullfile(root, "DATA_CONTRACT.md"); ...
             fullfile(root, "RELEASE_NOTES.md"); fullfile(root, "lsqvce_theory.pdf")];
    opts = matlab.addons.toolbox.ToolboxOptions(root, "lsqvce", ToolboxFiles=paths);
    opts.ToolboxName    = "lsqvce";
    opts.ToolboxVersion = string(ver{1});
    opts.Summary        = "Least-squares adjustment with variance component estimation";
    opts.ToolboxMatlabPath = root;
    if ~isfolder(fullfile(root, 'release')), mkdir(fullfile(root, 'release')); end
    opts.OutputFile = fullfile(root, 'release', "lsqvce_" + ver{1} + ".mltbx");
    matlab.addons.toolbox.packageToolbox(opts);
    fprintf('wrote %s\n', opts.OutputFile);
end
