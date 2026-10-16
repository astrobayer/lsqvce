function blocks = echoWarnings(capt)
%ECHOWARNINGS  Re-emit the MATLAB warnings of an inner evalc capture.
%   blocks = echoWarnings(capt)
%
%   WHY THIS EXISTS (WP-145, finding A1 of HANDOVER 12.77). A test that runs
%   code under EVALC to keep the console quiet also swallows every warning
%   that code raises: the warning is printed into the INNER capture and
%   discarded with it, so RUN_ALL_LSQVCE's outer capture -- the text its
%   warning inventory parses, and GATEVERDICT rule 7 judges -- never sees
%   it. Measured on main e483142: four of the fifteen examples raised
%   lsqvce:adjustVCEjoint:notConverged inside test_examples, and the gate
%   was green with an inventory of the probe alone. Every evalc in a suite
%   of dev/tests therefore hands its capture to this function, which prints
%   the warnings back out; test_ci_gate refuses an evalc that does not.
%
%   WHAT COUNTS AS A WARNING. Only the blocks MATLAB itself frames: in a
%   capture a warning is written as [\b Warning: <message> (Type "warning off
%   <id>" ...) ]\b, the frame being the two characters '[' char(8) ... ']'
%   char(8) (measured, also in a nested evalc, R2026a). A line that merely
%   STARTS with "Warning:" is text: demo_messaging replays its log sink on
%   screen, and copying those lines would inventory them as <unidentified>.
%   The block is re-emitted verbatim, frame and identifier hint included, on
%   a line of its own, so the outer parser reads it exactly as it reads a
%   warning raised at its own level.
%
%   INPUT
%     capt    1xN char/string   the text an evalc call returned
%
%   OUTPUT
%     blocks  1xW cell of char  the re-emitted blocks, in order ({} when the
%                               capture held no warning). Printed to the
%                               command window as a side effect; the return
%                               value is for the test of this function.
%
%   ERRORS
%     lsqvce:echoWarnings:badCapture   capt is not a char row or string scalar
%
%   EXAMPLE
%     out = evalc('[R, s] = demo_priors();');
%     echoWarnings(out);              % a warning of demo_priors reaches the gate
%
%   See also RUN_ALL_LSQVCE, GATEVERDICT, MUTEWARN
%
%   Created 2026-10-16 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

    narginchk(1, 1);
    if isstring(capt) && isscalar(capt), capt = char(capt); end
    if ~ischar(capt) || (~isempty(capt) && size(capt, 1) ~= 1)
        error('lsqvce:echoWarnings:badCapture', ...
              'capt must be the char row or string scalar an evalc call returned.');
    end
    bs = char(8);
    blocks = regexp(capt, ['\[' bs 'Warning:.*?\]' bs], 'match');
    for k = 1:numel(blocks)
        fprintf('\n%s\n', blocks{k});
    end
end
