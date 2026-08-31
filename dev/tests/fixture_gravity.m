function [cfg, L, data] = fixture_gravity()
%FIXTURE_GRAVITY  ENGINE-tier loaders over the demo_gravity_dataLoader data.
%
%   [cfg, L, data] = fixture_gravity()
%
%   The physically structured gravity fixture for the ENGINE-tier tests
%   (test_physOracle): since L10 (WP-74) the example returns the USER-tier
%   declaration only ([cfg, data], v3 contract) and this fixture builds the
%   Loaders from it -- the same makeLoader mapping LSQVCE.FIT performs
%   internally, written out once, where engine calls are legitimate.
%
%   INPUT   (none; examples/ must be on the path -- the callers add it)
%   OUTPUT
%     cfg   1x1 struct  the demo's configuration and planted truth
%     L     1x3 cell    Loaders, one per system (LSQVCE.ENGINE.MAKELOADER)
%     data  1x1 struct  the v3 data struct, passed through
%
%   ERRORS  none of its own; makeLoader construction errors pass through.
%
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- Level L10, WP-74.
    [cfg, data] = demo_gravity_dataLoader();
    L = cell(1, numel(data.systems));
    for k = 1:numel(data.systems)
        s = data.systems(k);
        L{k} = lsqvce.engine.makeLoader(s.buildDay, ...
                 'ObservationTypes', fieldnames(s.covariance).', ...
                 'Covariance',       s.covariance, ...
                 'Arcs',             s.arcs, ...
                 'LocalsPerArc',     s.localsPerArc, ...
                 'SystemId',         s.id, ...
                 'Name',             s.name, ...
                 'CommonParams',     data.commonParams, ...
                 'Days',             data.days);
    end
end
