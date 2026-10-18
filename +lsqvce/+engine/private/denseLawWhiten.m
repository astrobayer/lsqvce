function [A, B, l] = denseLawWhiten(C, A, B, l)
%DENSELAWWHITEN  Whiten a batch by the lower Cholesky factor of its dense law.
%   [A, B, l] = denseLawWhiten(C, A, B, l)
%
%   The ONE dense whitening step of the engine (WP-155, review 12.84 P1):
%   with C = L L' it returns L \ A, L \ B, L \ l -- forward substitutions on
%   the batch's columns, O(n^2 K). It used to form the explicit inverse
%   W = L \ eye(n) (LSQVCE.KERNEL.DECORRFACTOR) and multiply, O(n^3): two
%   thirds of a dense whitening's time at n = 2048 (0.037 s of the inverse
%   against 0.011 s of the factor), paid three times per arc and pass
%   (accumulation, closing solve, residual pass). The numbers are the same
%   to rounding (W A and L \ A agree to 1e-15). A C that Cholesky refuses
%   still takes decorrFactor's eigenvalue rescue, with its warning.
%
%   INPUT   C  n x n   the batch's dense covariance (symmetric)
%           A, B, l    the batch's rows (B may be sparse; the result is full,
%                      as W * B was)
%   OUTPUT  A, B, l    whitened
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-155.
    [L, fail] = chol(C, 'lower');                    % decorrFactor's own factor (lower triangle read)
    if fail == 0
        if ~isempty(A), A = L \ A; end
        if ~isempty(B), B = L \ B; end
        l = L \ l;
        return
    end
    W = lsqvce.kernel.decorrFactor(C);               % the rescue path (warns)
    A = W * A;  B = W * B;  l = W * l;
end
