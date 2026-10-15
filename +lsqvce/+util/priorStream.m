classdef priorStream
%LSQVCE.UTIL.PRIORSTREAM  The prior pseudo-streams' name: made and recognised in one place.
%
%   obsType = lsqvce.util.priorStream.of(groupName)
%   tf      = lsqvce.util.priorStream.is(obsTypes)
%
%   A local prior group of a system (data.systems(k).priors) enters the
%   adjustment as a PSEUDO-STREAM of that system: rows A = 0, B = L/sigma on
%   their day's local block, with their own law, estimated like a data
%   stream when the group says estimate "scale" (R5-D8 round 5). Its
%   observation type is 'prior_' followed by the group's name made a legal
%   field name. Every reader that treats those rows differently -- the
%   robust passes never re-weight them, snoop never rejects them, they get
%   no arc factor, before the group's 'from' pass they are inactive rather
%   than a data gap -- recognises the stream by that prefix, and a data
%   stream may not carry it (lsqvce:makeLoader:priorType).
%
%   WP-138 (code review of WP-135): the prefix was spelled out at nine
%   places -- startsWith(..., 'prior_') in the loader, the joint driver, the
%   robust residuals and the verdict, and fit parsed the stream key for
%   '|prior_'. It is defined HERE, once; the callers ask.
%
%   METHODS (static)
%     obsType = of(groupName)
%         INPUT   groupName  1xN char | 1x1 string  the group's name
%         OUTPUT  obsType    1xK char               'prior_<validName>'
%                 (matlab.lang.makeValidName: 'gravity field' -> 'gravityField')
%     tf = is(obsTypes)
%         INPUT   obsTypes   char | string | cellstr  observation type(s)
%         OUTPUT  tf         logical   1x1 for a char row, else the size of
%                                      obsTypes: true for a prior pseudo-stream
%
%   ACCURACY
%     Not applicable: string formatting and a prefix test.
%
%   ERRORS
%     none of its own
%
%   EXAMPLE
%     ot = lsqvce.util.priorStream.of('gravity field')   % 'prior_gravityField'
%     lsqvce.util.priorStream.is({ot, 'kbr'})             % [true false]
%
%   See also LSQVCE.UTIL.STREAMKEY, LSQVCE.ENGINE.MAKELOADER, LSQVCE.LAYOUT
%
%   Created 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   one definition of the prior pseudo-streams' prefix for the nine sites
%   that spelled it out. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.

    properties (Constant, Access = private)
        Prefix = 'prior_'      % the ONE spelling
    end

    methods (Static)
        function obsType = of(groupName)
        %OF  The observation type of a prior group's pseudo-stream.
            obsType = [lsqvce.util.priorStream.Prefix, ...
                       matlab.lang.makeValidName(char(groupName))];
        end

        function tf = is(obsTypes)
        %IS  True for the observation type of a prior pseudo-stream.
            tf = startsWith(obsTypes, lsqvce.util.priorStream.Prefix);
        end
    end
end
