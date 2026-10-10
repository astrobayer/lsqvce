function seedFor(tc, seed)
%SEEDFOR  Seed the global RNG for one test and restore it afterwards.
%   seedFor(tc, seed)
%
%   THE one way a test that draws from the GLOBAL stream seeds it (R7-45):
%   rng(seed) now, the previous generator state restored via tc.addTeardown,
%   so the inputs of the test do not depend on which suites ran before it,
%   and the suites after it do not depend on how many numbers it drew.
%   Called from setupOnce it scopes the whole suite (shared-fixture
%   teardown). Statistical tests use private RandStream objects instead;
%   this is for points whose inputs are arbitrary but must be repeatable.
%
%   INPUT
%     tc    1x1 matlab.unittest.TestCase   the running test case
%     seed  1x1 double                     nonnegative integer seed, as rng(seed)
%
%   OUTPUT  none (global RNG state altered until teardown)
%
%   ERRORS
%     Input:
%       lsqvce:seedFor:badSeed   seed is not a nonnegative integer scalar
%
%   EXAMPLE
%     seedFor(tc, 228);
%     raw = sign(randn(4, 10));               % repeatable in any suite order
%
%   See also MUTEWARN, RNG
%
%   Created 2026-10-10 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-117 (finding R7-45).

    narginchk(2, 2);
    if ~(isnumeric(seed) && isscalar(seed) && seed >= 0 && seed == round(seed))
        error('lsqvce:seedFor:badSeed', 'seed must be a nonnegative integer scalar.');
    end
    prev = rng;
    rng(seed);
    tc.addTeardown(@() rng(prev));
end
