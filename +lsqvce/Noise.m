classdef Noise
%LSQVCE.NOISE  The temporal law and level of one stream's noise.
%
%   n = lsqvce.Noise.white(sigma)            white noise
%   n = lsqvce.Noise.arma(ar, sigma)         a known whitening filter, MA = 1
%   n = lsqvce.Noise.arma(ar, sigma, MA=ma)
%   n = lsqvce.Noise.psd(f, S)               a measured one-sided spectrum, f from 0 Hz
%   n = lsqvce.Noise.acf(@rho, sigma)        correlation as a function of the lag
%   n = lsqvce.Noise.var(ar, Sigma)          vector AR: axes correlated in time
%   n = lsqvce.Noise.custom(@covFcn)         any dense covariance of an arc
%   n = lsqvce.Noise.<law>(..., MaxDenseArc=m)
%
%   One constructor per law replaces the struct('type', ...) specification of
%   a stream. The object states the LAW and its LEVEL only; what
%   Options.Estimate="noise" may change about it is a setting of the stream
%   (LSQVCE.SYSTEMDEFINITION/ADDSTREAM: Estimate, ArcFactors). The accuracy
%   of each epoch is not part of the law either: Build passes it with each
%   arc (LSQVCE.DAYDATA/ADDARC, Sigma=). sigma is always a STANDARD
%   DEVIATION, never a variance. Under the stream's Estimate="law" (the
%   default) only sigma is used, as the starting level, and the shape is
%   estimated from white noise.
%
%   LAWS
%     white(sigma)        sigma scalar, or a column of per-row sigmas (a fixed
%                         structure: its level can be estimated, its shape not)
%     arma(ar, sigma)     ar = [1 a1 ... ap], sigma the innovation standard
%                         deviation; MA = [1 b1 ... bq] (default 1), roots strictly
%                         inside the unit circle
%     psd(f, S)           f (Hz, one-sided, starting at 0) and S (unit^2/Hz),
%                         columns of one length; prepend f = 0 when the grid
%                         starts higher ([0; f], [S(1); S])
%     acf(@rho, sigma)    rho(dt, opt) the correlation at lag dt (seconds)
%     var(ar, Sigma)      ar d x d x p, Sigma d x d; no per-arc Sigma with it
%     custom(@covFcn)     C = covFcn(t, opt) the covariance of an arc at epochs t
%   NAME-VALUE (every law)
%     MaxDenseArc  positive integer   the longest arc factorised densely
%                                     (default 20 000 epochs; 3.2 GB of factor)
%   COST per arc of m rows: white O(m), arma O(m p), var O(m p d^2), psd O(m^3)
%   (Toeplitz) or O(m p) with Options.Noise.Whitener = "ar", acf / custom O(m^3).
%
%   INPUT   as the call forms above; every number finite
%   OUTPUT  n  1x1 lsqvce.Noise   properties Law (string), Parameters (struct),
%                                 MaxDenseArc ([] = the default)
%
%   ACCURACY
%     Not applicable: it stores a specification; the law's arithmetic is
%     LSQVCE.ENGINE.COVMODEL's. Each constructor equals the engines' covariance
%     specification field by field (test_taskDefinition).
%
%   ERRORS
%     MATLAB:validators:*                 a non-positive or non-finite sigma, a
%                                         handle that is not a function handle, a
%                                         non-integer MaxDenseArc
%     MATLAB:validation:UnableToConvert   a non-numeric argument
%     (the law's own consistency -- spectrum grid, MA invertibility -- is checked
%     by LSQVCE.ENGINE.COVMODEL when fit runs, under lsqvce:covModel:*)
%
%   EXAMPLE
%     lsqvce.Noise.white(1e-7)                       % let the data find the colour
%     lsqvce.Noise.arma([1 -1.6 0.9 -0.25], 1e-7)    % a known whitening filter
%     lsqvce.Noise.psd([0; 1e-3; 1e-2], [4e-14; 4e-14; 1e-14])
%
%   See also LSQVCE.SYSTEMDEFINITION, LSQVCE.TASKDEFINITION, LSQVCE.ENGINE.COVMODEL,
%            LSQVCE.ENGINE.CHECKPSD
%
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.

    properties (SetAccess = private)
        Law (1,1) string = "white"     % "white" | "arma" | "psd" | "acf" | "var" | "custom"
        Parameters (1,1) struct = struct()
        MaxDenseArc double = []
    end

    methods (Access = private)
        function obj = Noise(law, parameters, maxDenseArc)
            obj.Law = law;
            obj.Parameters = parameters;
            obj.MaxDenseArc = maxDenseArc;
        end
    end

    methods (Static)
        function n = white(sigma, opts)
            arguments
                sigma (:,1) double {mustBeNonempty, mustBePositive, mustBeFinite}
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("white", struct('sigma', sigma), opts.MaxDenseArc);
        end

        function n = arma(ar, sigma, opts)
            arguments
                ar (1,:) double {mustBeNonempty, mustBeFinite}
                sigma (1,1) double {mustBePositive, mustBeFinite}
                opts.MA (1,:) double {mustBeNonempty, mustBeFinite} = 1
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("arma", struct('ar', ar, 'ma', opts.MA, 'sigma', sigma), opts.MaxDenseArc);
        end

        function n = psd(f, S, opts)
            arguments
                f (:,1) double {mustBeNonempty, mustBeNonnegative, mustBeFinite}
                S (:,1) double {mustBeNonempty, mustBeNonnegative, mustBeFinite}
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("psd", struct('f', f, 'S', S), opts.MaxDenseArc);
        end

        function n = acf(rho, sigma, opts)
            arguments
                rho (1,1) function_handle
                sigma (1,1) double {mustBePositive, mustBeFinite}
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("acf", struct('acf', rho, 'sigma', sigma), opts.MaxDenseArc);
        end

        function n = var(ar, Sigma, opts)
            arguments
                ar (:,:,:) double {mustBeNonempty, mustBeFinite}
                Sigma (:,:) double {mustBeNonempty, mustBeFinite}
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("var", struct('ar', ar, 'sigma', Sigma), opts.MaxDenseArc);
        end

        function n = custom(covFcn, opts)
            arguments
                covFcn (1,1) function_handle
                opts.MaxDenseArc double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
            end
            n = lsqvce.Noise("custom", struct('fcn', covFcn), opts.MaxDenseArc);
        end
    end

    methods (Hidden)
        function s = spec(obj)
            % SPEC  The engines' covariance specification of this law (type + fields).
            types = struct('white', 'diag', 'arma', 'arma', 'psd', 'psd', 'acf', 'acf', ...
                           'var', 'var', 'custom', 'custom');
            s = struct('type', types.(char(obj.Law)));
            f = fieldnames(obj.Parameters);
            for k = 1:numel(f), s.(f{k}) = obj.Parameters.(f{k}); end
            if ~isempty(obj.MaxDenseArc), s.maxDenseArc = obj.MaxDenseArc; end
        end
    end
end
