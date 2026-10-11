function [AtilBar, Bbar, ebar] = jointArcParts(A, B, l, W, x, y, G)
%JOINTARCPARTS  Decorrelated joint-arc pieces for the per-lag estimators.
%   [AtilBar, Bbar, ebar] = jointArcParts(A, B, l, W, x, y, G)
%
%   The ONE shared implementation (rule 3) behind the two joint spectral
%   sweeps (adjustVCEjoint's oneArc; until M4 4c also the separate mode's) that
%   makes both consistent with Ellmer's step 7 (the FULL Table-6.2 parameter
%   vector) and with the joint projector (WP-33, fixes finding F-03):
%
%     AtilBar = W^-T (A - B G)      the REDUCED common design, decorrelated.
%                                   Its hat against the reduced normal
%                                   N = Ncc_red is exactly the common part of
%                                   the JOINT hat (block identity, see
%                                   LSQVCE.ENGINE.BLOCKCOFACTOR), so the shared probe
%                                   U^-1 Zbar built from chol(N) is finally
%                                   probing the operator it factors.
%     Bbar    = W^-T B              the decorrelated local design; the local
%                                   hat B Nll^-1 B' is added EXACTLY (nl
%                                   columns) in LSQVCE.KERNEL.REDUNDANCYLAGS via its
%                                   'LocalBbar'/'LocalLl' inputs.
%     ebar    = W^-T (l - A x - B y)  the residual of the FULL joint fit --
%                                   NOT l - A x: omitting B y leaves the
%                                   arc-parameter signal in the "noise"
%                                   (measured 3.8x Omega_0 inflation,
%                                   v10review.py T3).
%
%   INPUT
%     A   m x nc double   the arc's RAW common design (as the loader hands it)
%     B   m x nl double   the arc's RAW local design; may have nl = 0 columns
%     l   m x 1  double   the arc's raw observations
%     W   m x m  double   the arc's whitening factor, upper triangular with
%                         Sigma = W'*W (as required by LSQVCE.KERNEL.REDUNDANCYLAGS)
%     x   nc x 1 double   the current COMMON solution (pass sol.x(1:nc))
%     y   nl x 1 double   this SYSTEM's local solution, recovered once per
%                         pass as Nll_r \ (nl_r - Nlc*x) -- exactly as
%                         LSQVCE.ENGINE.REDUNDANCY does
%     G   nl x nc double  this SYSTEM's Schur factor inv(Nll) * Nlc (the same
%                         G as LSQVCE.ENGINE.BLOCKCOFACTOR)
%
%   OUTPUT
%     AtilBar  m x nc double   W^-T (A - B G)
%     Bbar     m x nl double   W^-T B  (m x 0 when nl = 0)
%     ebar     m x 1  double   W^-T (l - A x - B y)
%
%   nl = 0 (no local parameters): B has zero columns, y and G are empty --
%   then AtilBar = W^-T A and ebar = W^-T (l - A x), i.e. EXACTLY the
%   pre-WP-33 quantities, so the nl = 0 path is unchanged by construction.
%
%   ERRORS
%     lsqvce:jointArcParts:sizeMismatch  A, B, l or x/y/G disagree in size.
%
%   EXAMPLE
%     y = st.Nll_r \ (st.nl_r - st.Nlc * x);      % once per system and pass
%     G = st.Nll_r \ st.Nlc;
%     [At, Bb, eb] = jointArcParts(a.A, a.B, a.l(:), W, x, y, G);
%     r = lsqvce.kernel.redundancyLags(At, W, UiZ, eb, nLag, ...
%                               'LocalBbar', Bb, 'LocalLl', Ll);
%
%   See also LSQVCE.KERNEL.REDUNDANCYLAGS, LSQVCE.ENGINE.BLOCKCOFACTOR, LSQVCE.ENGINE.REDUNDANCY,
%            LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): the
%   retired leverage references point at lsqvce.engine.BlockCofactor. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    m  = size(A, 1);
    nc = size(A, 2);
    nl = size(B, 2);
    if size(B,1) ~= m || numel(l) ~= m
        error('lsqvce:jointArcParts:sizeMismatch', ...
              'A (%dx%d), B (%dx%d) and l (%d) must share their row count.', ...
              m, nc, size(B,1), nl, numel(l));
    end
    if numel(x) < nc
        error('lsqvce:jointArcParts:sizeMismatch', ...
              'x has %d elements but A has %d columns.', numel(x), nc);
    end
    x = x(:);
    l = l(:);

    Abar = W.' \ A;
    if nl == 0
        Bbar    = zeros(m, 0);
        AtilBar = Abar;
        ebar    = W.' \ (l - A * x(1:nc));
        return
    end
    if numel(y) ~= nl || size(G,1) ~= nl || size(G,2) ~= nc
        error('lsqvce:jointArcParts:sizeMismatch', ...
              'y (%d) / G (%dx%d) do not match nl = %d, nc = %d.', ...
              numel(y), size(G,1), size(G,2), nl, nc);
    end
    Bbar    = W.' \ B;
    AtilBar = Abar - Bbar * G;                 % W^-T (A - B G)
    ebar    = W.' \ (l - A * x(1:nc) - B * y(:));
end
