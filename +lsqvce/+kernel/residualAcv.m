function [acv, nSamp] = residualAcv(res, satId, obsType, nAxes, maxLag, dayFcn)
%LSQVCE.KERNEL.RESIDUALACV  Temporal autocovariance of de-scaled residuals.
%   [acv, nSamp] = lsqvce.kernel.residualAcv(res, satId, obsType, nAxes, maxLag, dayFcn)
%
%   Estimates the TEMPORAL autocovariance of one satellite/obsType stream after
%   removing the FIXED per-epoch factor D that the whitener already applied. It
%   pools the ACV over the arcs and over the axes (the axes share one temporal
%   law R_t, the separable model C = D (R_t kron I_d) D). Pooling never crosses
%   satellites: pass a single satId.
%
%   De-scaling by D:
%     nAxes == 1 (e.g. rr)   divide each epoch residual by dayData.sigma(rows)
%     nAxes  > 1 (e.g. pos)  left-multiply each epoch residual by
%                            chol(dayData.Cep(epoch))^-1 (lower)
%
%   INPUT
%     res      1x1   struct   Collected residuals (see LSQVCE.ENGINE.COLLECTRESIDUALS).
%                            Uses res.sys(:).sysID and res.sys(:).obs(:).
%                            {obsType, day, arc, v}. Required.
%     satId    1x1   char|string|numeric   Satellite to select. Required.
%     obsType  1xk   char|string   Observation type, e.g. 'pos'. Required.
%     nAxes    1x1   integer >=1   Number of axes per epoch (rr=1, pos=3).
%                            Required.
%     maxLag   1x1   integer >=1   Highest ACV lag to accumulate. Required.
%     dayFcn   1x1   function_handle   dayFcn(satId, day) -> day struct exposing
%                            .arc(iArc).rows / .epochs and .sigma (nAxes==1) or
%                            .Cep (nEpochs x nAxes x nAxes). Must already know
%                            obsType (bake it into the closure). Required.
%
%   OUTPUT
%     acv      (maxLag+1)x1   double   Biased (1/n) autocovariance r(0..maxLag);
%                            biased so it stays positive semidefinite and
%                            LSQVCE.KERNEL.LEVINSONDURBIN never trips on it.
%     nSamp    1x1           double    Axis-sample count pooled into acv; feed to
%                            LSQVCE.ENGINE.SELECTARACV for its AICc/BIC/PACF reporting.
%
%   EXAMPLE
%     dayFcn = @(s,d) myLoader(cfg, s, d, 'pos', 'day');   % your day builder
%     [r, n] = lsqvce.kernel.residualAcv(res, 'C', 'pos', 3, 128, dayFcn);
%     [a, s2, p] = lsqvce.engine.selectArAcv(r, n);
%
%   ERRORS
%     lsqvce:residualAcv:noObs      the stream carries no observations
%     lsqvce:residualAcv:dayFcn     the day accessor is not callable
%
%   See also LSQVCE.ENGINE.SELECTARACV, LSQVCE.ENGINE.COLLECTRESIDUALS,
%            LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.SEPARABLEWHITEN
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
    narginchk(6, 6);
    validateattributes(nAxes,  {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nAxes', 4);
    validateattributes(maxLag, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'maxLag', 5);
    if ~isa(dayFcn, 'function_handle')
        error('lsqvce:residualAcv:dayFcn', 'dayFcn must be a function handle.');
    end

    acvNumerator = zeros(maxLag + 1, 1);
    nSamp        = 0;
    dayCache     = containers.Map('KeyType','char', 'ValueType','any');

    for iSys = 1:numel(res.sys)
        if ~isequal(res.sys(iSys).sysID, satId), continue, end   % this satellite
        arcObs = res.sys(iSys).obs;
        for iArc = 1:numel(arcObs)
            if ~strcmp(arcObs(iArc).obsType, obsType), continue, end

            key = sprintf('%s|%d', string(satId), arcObs(iArc).day);
            if ~dayCache.isKey(key)
                dayCache(key) = dayFcn(satId, arcObs(iArc).day);
            end
            dayData = dayCache(key);
            arc     = dayData.arc(arcObs(iArc).arc);

            v      = arcObs(iArc).v(:);              % raw residuals, epoch-major
            nE     = numel(v) / nAxes;
            byAxis = reshape(v, nAxes, nE);          % axis x epoch

            % ---- de-scale by the FIXED D --------------------------------
            if nAxes == 1
                epochStd = dayData.sigma(arc.rows);
                byAxis   = byAxis ./ epochStd(:).';
            else
                Cep = dayData.Cep(arc.epochs, :, :);
                for iE = 1:nE
                    Ck = reshape(Cep(iE,:,:), nAxes, nAxes);
                    Ck = (Ck + Ck.')/2;
                    byAxis(:,iE) = chol(Ck,'lower') \ byAxis(:,iE);
                end
            end

            % ---- pool ACV over the axes ---------------------------------
            for iAxis = 1:nAxes
                s = byAxis(iAxis,:).';  s = s - mean(s);
                m = numel(s);  L = min(maxLag, m - 1);
                for lag = 0:L
                    acvNumerator(lag+1) = acvNumerator(lag+1) ...
                        + s(1:m-lag).' * s(1+lag:m);
                end
                nSamp = nSamp + m;
            end
        end
    end
    if nSamp == 0
        error('lsqvce:residualAcv:noObs', ...
              'No %s residuals for satellite %s.', obsType, string(satId));
    end
    acv = acvNumerator / nSamp;                      % biased (1/n): stays PSD
end
