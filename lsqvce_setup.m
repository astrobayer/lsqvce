function lsqvce_setup()
%LSQVCE_SETUP  Add lsqvce to the MATLAB path and register its Help browser docs.
%
%   Run this once from the toolbox root. It:
%     * adds the toolbox root (so the +lsqvce package is visible),
%     * adds tests/, examples/ and template/,
%     * rebuilds the Help browser index so `doc lsqvce` works.
%
%   After running:
%       doc lsqvce                 % the toolbox landing page
%       doc lsqvce.diagnose        % any individual function
%       help lsqvce.snoop          % command-window help
%       run_all_lsqvce              % the full test suite (THE gate)
%       demo_workflow               % the interpretation walk-through
%       edit lsqvce_template        % start your own case
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

root = fileparts(mfilename('fullpath'));

addpath(root);
addpath(fullfile(root, 'tests'));
addpath(fullfile(root, 'examples'));
addpath(fullfile(root, 'template'));
% One flat +lsqvce package: adding the root makes every lsqvce.* function
% visible. tests/, examples/ and template/ are folder names, not a namespace
% tier.
%
% DO NOT add a folder here without creating it: addpath warns
% (MATLAB:addpath:nameNonexistentOrNotADirectory) on a missing target, and
% three dead entries (tests/accum, examples/accum, examples/gravity, all
% removed 2026-07-18) sat here unnoticed until 2026-07-19 because nothing
% checked. audit_tree.py check 12 now does.

fprintf('lsqvce: added to path.\n');

% The Help browser reads info.xml from a folder ON THE PATH, and helptoc.xml
% from the location info.xml points at. Both live in the toolbox root.
if exist(fullfile(root, 'info.xml'), 'file')
    try
        builddocsearchdb(root);   % root holds helptoc.xml + info.xml; the
        fprintf('lsqvce: help search database built.\n');
    catch ME
        fprintf(['lsqvce: could not build the doc search DB (%s).\n' ...
                 '         The HTML pages still work: open html/index.html, ' ...
                 'or use `help lsqvce.foo`.\n'], ME.message);
    end
end

fprintf('\n  doc lsqvce              the landing page (User Manual, guides, API)\n');
fprintf('  demo_workflow            build your OWN case, from scratch\n');
fprintf('  demo_lsqvce             six failure modes\n');
fprintf('  lsqvce_template          start your OWN case (template/)\n');
fprintf('  run_all_lsqvce           the unit-test suite (THE gate)\n\n');
fprintf('  All docs are HTML under doc lsqvce: User/Reference Manual,\n');
fprintf('  Interpretation Guide, and "Choosing an engine".\n\n');

if ~savepath()
    fprintf('lsqvce: path saved permanently.\n');
else
    fprintf(['lsqvce: could not save the path (no write permission). ' ...
             'Re-run lsqvce_setup each session, or add it to startup.m.\n']);
end
end
