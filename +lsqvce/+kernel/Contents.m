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
%   residualAcv          - autocovariance of residuals
%   separableWhiten      - separable (temporal x axis) whitening
%   toeplitzCofactor     - Toeplitz cofactor construction
%   varACV               - VAR autocovariance
%   varWhiten            - VAR whitening
%   whitenessCorr        - whiteness (Ljung-Box) statistic
