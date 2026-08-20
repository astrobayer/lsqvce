function n = declaredPoints(file)
%DECLAREDPOINTS  Count the test points DECLARED in a function-based suite file.
%   n = declaredPoints(file)
%
%   Counts lines matching the test-point signature
%   `function [out =] test<Name>(...)`, excluding the suite entry point
%   `function tests = ...` and anything inside a %{ ... %} block comment.
%   RUN_ALL_LSQVCE compares this against the EXECUTED count per suite
%   (WP-47, adopted from the external geoMap handover 2.6.4): a suite that
%   silently drops points -- a load-time filter, a shadowed local, a
%   block-commented test -- is indistinguishable from green by the pass
%   count alone. The Python audit (check 16) counts the same definition
%   independently; the two counters deliberately share NO code (the
%   second gate reads something the first does not, HANDOVER 6.2.3).
%
%   INPUT
%     file  1xK char/string  absolute path of a tests/test_*.m file.
%
%   OUTPUT
%     n     1x1 double       number of declared test points.
%
%   ERRORS
%     File access:
%       lsqvce:declaredPoints:noFile   file does not exist or is unreadable
%
%   EXAMPLE
%     n = declaredPoints(fullfile(root, 'tests', 'test_util.m'));
%
%   See also RUN_ALL_LSQVCE, GATEVERDICT
%
%   Created 2026-08-20 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    narginchk(1, 1);
    file = char(file);
    if exist(file, 'file') ~= 2
        error('lsqvce:declaredPoints:noFile', 'no such file: %s', file);
    end
    txt = fileread(file);
    lines = strsplit(txt, '\n', 'CollapseDelimiters', false);

    n = 0;
    inBlock = false;
    for k = 1:numel(lines)
        s = strtrim(lines{k});
        % Block comments: a %{ / %} line must be ALONE on its line to open/
        % close (MATLAB semantics); nesting is tracked as a flag because the
        % suites do not nest them.
        if strcmp(s, '%{'), inBlock = true;  continue; end
        if strcmp(s, '%}'), inBlock = false; continue; end
        if inBlock, continue; end
        if ~isempty(regexp(s, '^function\s+tests\s*=', 'once'))
            continue                                   % the suite entry point
        end
        if ~isempty(regexp(s, ...
                '^function\s+(?:\w+\s*=\s*)?test\w+\s*\(', 'once'))
            n = n + 1;
        end
    end
end
