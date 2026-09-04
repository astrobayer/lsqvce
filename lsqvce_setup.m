function lsqvce_setup(action)
%LSQVCE_SETUP  Install lsqvce (path + Help browser docs), check it, or uninstall it.
%
%   lsqvce_setup                install: add the toolbox to the path, register the
%                               Help browser documentation, save the path, verify
%   lsqvce_setup('install')     the same
%   lsqvce_setup('check')       report what is installed (path entries, doc index,
%                               version, whether `doc lsqvce` and `lsqvce.fit` resolve)
%   lsqvce_setup('uninstall')   remove the toolbox from the path, save the path, remove
%                               the Help browser search database it built
%
%   WHAT INSTALL DOES (and nothing else)
%     * addpath(root)                 the lsqvce.* packages
%     * addpath(root/examples)        the demos (demo_fit, demo_gravfire, ...)
%     * builddocsearchdb(root)        the Help browser index: `doc lsqvce` opens the
%                                     Manual (html/lsqvce_manual.html); info.xml and
%                                     helptoc.xml live in the root
%     * savepath                      so the next MATLAB session has it (a warning, not
%                                     an error, when pathdef.m is not writable)
%   Not on the path by design: dev/tests (the suites manage their path themselves),
%   template/ (copied into your project, not walked), dev/ (development records).
%
%   WHAT UNINSTALL DOES
%     rmpath of exactly the two entries install added; savepath; deletes the
%     helpsearch-v4 database folder builddocsearchdb created under the root. Your
%     own copies of template/gravfire_data.m are yours and stay where they are.
%
%   INPUT   action  char/string  'install' (default) | 'check' | 'uninstall'
%   OUTPUT  none (prints a report)
%   ERRORS  lsqvce:setup:badAction   an unknown action word
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, the Manual (doc lsqvce)
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8, 2026-07-11.
%   Rewritten 2026-09-04 (v5.2 docs level) with assistance from Claude (Anthropic),
%   model Claude Fable 5.1: install / check / uninstall as one entry point.
arguments
    action (1,1) string = "install"
end
    root = fileparts(mfilename('fullpath'));
    entries = {root, fullfile(root, 'examples')};
    switch lower(action)
        case "install"
            for k = 1:numel(entries)
                if ~isfolder(entries{k})
                    error('lsqvce:setup:badAction', 'expected folder missing: %s', entries{k});
                end
                addpath(entries{k});
            end
            fprintf('lsqvce: added to the path (root, examples).\n');
            if isfile(fullfile(root, 'info.xml'))
                try
                    builddocsearchdb(root);
                    fprintf('lsqvce: Help browser search database built.\n');
                catch ME
                    fprintf(['lsqvce: could not build the doc search DB (%s).\n' ...
                             '        The HTML pages still work: open html/lsqvce_manual.html, ' ...
                             'or use `help lsqvce.fit`.\n'], ME.message);
                end
            end
            if savepath() ~= 0
                warning('lsqvce:setup:savepath', ...
                        'the path could not be saved (pathdef.m not writable); this session has it.');
            end
            lsqvce_setup("check");
        case "check"
            v = lsqvce.util.toolboxVersion();
            fprintf('\nlsqvce %s at %s\n', v, root);
            for k = 1:numel(entries)
                onPath = any(strcmpi(strsplit(path, pathsep), entries{k}));
                fprintf('  path entry %-52s %s\n', entries{k}, tern(onPath, 'ok', 'MISSING'));
            end
            fprintf('  lsqvce.fit resolves                                       %s\n', tern(exist('lsqvce.fit', 'file') == 2, 'ok', 'MISSING'));
            fprintf('  Help browser index (helpsearch-v4)                        %s\n', tern(isfolder(fullfile(root, 'helpsearch-v4')), 'ok', 'not built'));
            fprintf('  manual html/lsqvce_manual.html                            %s\n', tern(isfile(fullfile(root, 'html', 'lsqvce_manual.html')), 'ok', 'MISSING'));
            fprintf('\n  doc lsqvce                    the Manual (the single authority)\n');
            fprintf('  help lsqvce.fit               command-window help of the one verb\n');
            fprintf('  demo_workflow                 build your own case, step by step\n');
            fprintf('  template/gravfire_data.m      copy into your project to start a GRACE case\n');
            fprintf('  lsqvce_setup(''uninstall'')     remove it again\n\n');
        case "uninstall"
            for k = 1:numel(entries)
                if any(strcmpi(strsplit(path, pathsep), entries{k}))
                    rmpath(entries{k});
                end
            end
            db = fullfile(root, 'helpsearch-v4');
            if isfolder(db)
                try, rmdir(db, 's'); catch, end
            end
            if savepath() ~= 0
                warning('lsqvce:setup:savepath', 'the path could not be saved; this session is clean.');
            end
            fprintf('lsqvce: removed from the path and the Help browser index; the folder %s is untouched.\n', root);
        otherwise
            error('lsqvce:setup:badAction', 'action must be ''install'', ''check'' or ''uninstall'' (got ''%s'').', action);
    end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
