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
%   hatACV               - ACV of every bin of a piecewise-linear spectrum
%   levinsonDurbin       - Levinson-Durbin recursion
%   psd2acv              - PSD -> autocovariance
%   lawInverse           - explicit inverse of an arc's covariance law
%   redundancyLags       - redundancy per lag without the dense W^-1
%   truncationDeficit    - the redundancy a robust sample loses to its truncation
%   wSnoop               - Baarda's w per epoch with iterative data snooping
%   separableWhiten      - separable (temporal x axis) whitening
%   stat                 - distributions and test thresholds, base MATLAB only
%   toeplitzCofactor     - Toeplitz cofactor construction
%   varACV               - VAR autocovariance
%   varWhiten            - VAR whitening
%   whitenessCorr        - whiteness (Ljung-Box) statistic
%   whittle            - multichannel Levinson (Whittle) VAR fit
%   robustWeights      - IRLS weight factors from standardised residuals
