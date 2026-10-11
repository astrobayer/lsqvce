% lsqvce.kernel -- the KERNEL tier of lsqvce (R4-D1, WP-64a).
%
% Users reach these through lsqvce.fit / lsqvce.Options / lsqvce.Results;
% experts call them directly. Roles from docs_src/tiers.csv, the ONE tier
% authority (audit check 18).
%
%   ar2arma              - AR -> ARMA conversion
%   armaACV              - ARMA autocovariance
%   armaWhiten           - ARMA whitening filter (batched)
%   crossCorr            - cross-correlation of residual streams
%   dctI                 - DCT-I transform
%   decorrFactor         - decorrelation factor of a dense covariance
%   levinsonDurbin       - Levinson-Durbin recursion
%   psd2acv              - PSD -> autocovariance
%   redundancyLags       - redundancy per lag without the dense W^-1
%   separableWhiten      - separable (temporal x axis) whitening
%   stat                 - distributions and test thresholds, base MATLAB only (from engine, WP-86)
%   toeplitzCofactor     - Toeplitz cofactor construction
%   varACV               - VAR autocovariance
%   varWhiten            - VAR whitening
%   whitenessCorr        - whiteness (Ljung-Box) statistic
%   whittle            - multichannel Levinson (Whittle) VAR fit (moved from engine, R5 M2)
%   robustWeights      - IRLS weight factors from standardised residuals (R5 M2)
