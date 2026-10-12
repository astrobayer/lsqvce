% lsqvce.util -- shared utilities (KERNEL tier; R4-D1, WP-64a).
%
%   callUser           - call a user function: outputs counted once, errors wrapped
%   getfielddef        - struct field with default
%   toolboxVersion     - the version string from Contents.m (one source)
%   onPass             - fire the per-pass inspection hook (Options.OnPass)
%   rowWeightKey       - the one key convention for per-arc robust weights
%   streamKey          - the one stream key 'satKey|obsType'
%   arcKey             - the one arc key of the joint estimator
%   log                - toolbox-wide logging state
%   mergeStruct        - struct merge (later wins)
%   verbosity          - verbosity level resolution and formatting
