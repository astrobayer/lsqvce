function [cfg, L, data] = fixture_grace()
%FIXTURE_GRACE  ENGINE-tier loaders over the demo_grace_dataLoader data.
%
%   [cfg, L, data] = fixture_grace()
%
%   The synthetic GRACE-FO fixture for the ENGINE-tier tests (test_adjustVCE):
%   since L10 (WP-74) the example returns the USER-tier declaration only
%   ([cfg, data], v3 contract) and this fixture builds the Loaders from it --
%   the same makeLoader mapping LSQVCE.FIT performs internally, written out
%   once, where engine calls are legitimate.
%
%   INPUT   (none; examples/ must be on the path -- the callers add it)
%   OUTPUT
%     cfg   1x1 struct  the demo's configuration and planted truth
%     L     1x2 cell    Loaders: the merged C/D pair and independent E
%     data  1x1 struct  the v3 data struct, passed through
%
%   ERRORS  none of its own; makeLoader construction errors pass through.
%
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): 'Arcs' is passed
%   only when the declaration carries a schedule. With assistance from
%   Claude (Anthropic), model Claude Fable 5.1.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- Level L10, WP-74.
    [cfg, data] = demo_grace_dataLoader();
    L = cell(1, numel(data.systems));
    for k = 1:numel(data.systems)
        s = data.systems(k);
        nv = {};                                           % v6: `arcs` is optional
        if isfield(s, 'arcs') && ~isempty(s.arcs), nv = {'Arcs', s.arcs}; end
        L{k} = lsqvce.engine.makeLoader(s.buildDay, nv{:}, ...
                 'ObservationTypes',  fieldnames(s.covariance).', ...
                 'Covariance',        s.covariance, ...
                 'SystemId',          s.id, ...
                 'Name',              s.name, ...
                 'LocalParams',       s.localParams, ...
                 'Scoped',            true, ...
                 'ProcessingOptions', s.processingOptions, ...
                 'CommonParams',      data.commonParams, ...
                 'Days',              data.days);
    end
end
