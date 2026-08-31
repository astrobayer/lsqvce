function f = makeKeyFcn(g)
%MAKEKEYFCN  Resolve a GroupBy spec into a key function f(sat,obsType,arc,day).
%   f = makeKeyFcn(g)
%
%   The ONE key-preset definition (WP-38, fixes F-08; rule 3) behind
%   LSQVCE.ENGINE.REDUNDANCY and LSQVCE.ENGINE.STREAMCONTRIBUTION -- previously two
%   verbatim copies that could drift apart, silently breaking the
%   cross-path group identity the streamed Hutchinson denominators depend
%   on (the 'arc' preset is also the keyed-arc identity of
%   LSQVCE.ENGINE.ADJUSTVCEJOINT, WP-34). ADJUSTVCE.vceGroupKeyFcn is a
%   DIFFERENT, class-scoped scheme and deliberately not unified.
%
%   INPUT
%     g   char | 1x1 function_handle   preset name ('arc' | 'obsType' |
%          'system' | 'day') or a custom @(sat,obsType,arc,day) -> char
%          handle (passed through unchanged)
%
%   OUTPUT
%     f   1x1 function_handle   @(sat,obsType,arc,day) -> group key char
%
%   ERRORS
%     lsqvce:makeKeyFcn:badGroupBy   g is neither a handle nor a known
%                                    preset name.
%
%   EXAMPLE
%     f = makeKeyFcn('arc');
%     key = f(71, 'pos', 3, 12);     % 's71_pos_d12_a3'
%
%   See also LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.STREAMCONTRIBUTION,
%            LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    if isa(g, 'function_handle')
        f = g;  return
    end
    switch lower(g)
        case 'arc',     f = @(s,ot,a,d) sprintf('s%s_%s_d%s_a%s', ...
                                    string(s), string(ot), string(d), string(a));
        case 'obstype', f = @(s,ot,a,d) sprintf('%s', string(ot));
        case 'system',  f = @(s,ot,a,d) sprintf('s%s', string(s));
        case 'day',     f = @(s,ot,a,d) sprintf('d%s', string(d));
        otherwise
            error('lsqvce:makeKeyFcn:badGroupBy', ...
                ['GroupBy must be ''arc'', ''obsType'', ''system'', ''day'' ' ...
                 'or a function handle (got ''%s'').'], g);
    end
end
