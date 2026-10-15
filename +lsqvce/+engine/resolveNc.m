function [nc, info] = resolveNc(systems)
%LSQVCE.ENGINE.RESOLVENC  The common-parameter count, from the systems themselves.
%
%   nc          = LSQVCE.ENGINE.RESOLVENC(systems)
%   [nc, info]  = LSQVCE.ENGINE.RESOLVENC(systems)
%
%   THE PARAMETRISATION IS THE SINGLE SOURCE OF TRUTH. nc used to be a leading
%   argument on LSQVCE.ENGINE.ADJUSTVCEJOINT and LSQVCE.ENGINE.ADJUSTBLOCK, which meant the same
%   number was stated twice -- once as the argument and once as the number of
%   rows in the 'commonParams' table -- free to disagree. The demos showed the
%   failure mode: `cfg.nc = 6` on one line and `cfg.common` with six rows on
%   another. This function removes the duplication by deriving nc from the
%   systems, so it can be stated only ONCE, where the system is defined.
%
%   WHERE nc COMES FROM, by system form
%     Loader object (LSQVCE.ENGINE.MAKELOADER)   L.nc, itself derived from the rows of
%                                         the 'commonParams' parametrisation.
%     Raw spec struct                     the .nc field, declared beside the .nl
%                                         the spec already carries.
%   nc is SHARED across systems by definition -- it is the common block every
%   system contributes to -- so all systems must agree, and disagreement is an
%   error rather than something silently resolved.
%
%   INPUT
%     systems  1xM cell | 1xM struct   loader objects, raw spec structs, or a
%                                      mix. See LSQVCE.ENGINE.ADJUSTBLOCK for the forms.
%
%   OUTPUT
%     nc    1x1 double   the common-parameter count
%     info  1x1 struct   .source   1xM cell, where each system's nc came from
%                        ('loader' | 'spec')
%                        .values   1xM double, each system's claim
%
%   ACCURACY
%     Not applicable: an integer count read from the systems and compared
%     exactly, no floating-point result (test_errorDoctrine
%     (testResolveNcContract)). WP-139.
%
%   ERRORS
%     lsqvce:resolveNc:empty        no systems given
%     lsqvce:resolveNc:missing      a system declares no nc at all
%     lsqvce:resolveNc:conflict     systems disagree
%
%   EXAMPLE
%     L  = lsqvce.engine.makeLoader(bd, 'commonParams', tbl, 'localParams', loc, ...);
%     nc = lsqvce.engine.resolveNc({L});          % = size(tbl,1), never restated
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    if isempty(systems)
        error("lsqvce:resolveNc:empty", ...
              "No systems were given, so the number of common parameters " + ...
              "cannot be determined. Pass at least one loader or spec.");
    end

    if isstruct(systems)
        systems = num2cell(systems);
    end
    if ~iscell(systems)
        systems = {systems};
    end

    nSys   = numel(systems);
    values = nan(1, nSys);
    source = cell(1, nSys);

    for k = 1:nSys
        s = systems{k};
        [values(k), source{k}] = ncOfOne(s, k);
    end

    good = values(isfinite(values));
    if isempty(good)
        error("lsqvce:resolveNc:missing", ...
              "None of the %d system(s) declares a common-parameter count. " + ...
              "Build them with LSQVCE.ENGINE.MAKELOADER using 'commonParams', or " + ...
              "add an .nc field to each raw spec struct.", nSys);
    end
    if any(good ~= good(1))
        error("lsqvce:resolveNc:conflict", ...
              "The systems disagree about the number of common parameters " + ...
              "(%s). nc is the SHARED block every system contributes to, so " + ...
              "it must be identical everywhere. Check that every loader was " + ...
              "built with the same 'commonParams' table.", ...
              join(string(good), ", "));
    end

    nc   = double(good(1));
    info = struct("source", {source}, "values", values);
end

%--------------------------------------------------------------------------
function [nc, src] = ncOfOne(s, k)
% Extract one system's declared nc, or NaN when it declares none.
    nc  = NaN;
    src = "none";
    if isstruct(s) && isfield(s, "isLoader") && isequal(s.isLoader, true)
        % A Loader from makeLoader: nc was derived from its commonParams rows.
        if isfield(s, "nc") && ~isempty(s.nc) && s.nc > 0
            nc  = s.nc;
            src = "loader";
        else
            error("lsqvce:resolveNc:missing", ...
                  "System %d is a Loader built WITHOUT a 'commonParams' " + ...
                  "parametrisation, so it cannot say how many common " + ...
                  "parameters there are. Rebuild it with 'commonParams'.", k);
        end
    elseif isstruct(s) && isfield(s, "nc") && ~isempty(s.nc)
        % A raw spec struct declaring nc beside the nl it already carries.
        nc  = s.nc;
        src = "spec";
    end
end
