% lsqvce.util -- shared utilities (KERNEL tier; R4-D1, WP-64a).
%
%   callUser           - call a user function: outputs counted once, errors wrapped (R5 M1)
%   getfielddef        - struct field with default
%   toolboxVersion     - the version string from Contents.m (one source) (v5.2)
%   onPass             - fire the per-pass inspection hook (Options.OnPass) (R5 M5)
%   rowWeightKey       - the one key convention for per-arc robust weights (R5 M2)
%   streamKey          - the one stream key 'satKey|obsType' (R7-43)
%   arcKey             - the one arc key of the joint estimator (R7-43)
%   log                - toolbox-wide logging state
%   mergeStruct        - struct merge (later wins)
%   verbosity          - verbosity level resolution and formatting
