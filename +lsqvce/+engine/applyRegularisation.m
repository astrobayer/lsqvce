function [Nr, nr, dg] = applyRegularisation(Nll, nl, dg, opt, mObs, ltl)
%LSQVCE.ENGINE.APPLYREGULARISATION  Regularise a local block before pre-elimination.
%
%   [Nr, nr, dg] = LSQVCE.ENGINE.APPLYREGULARISATION(Nll, nl, dg, opt, mObs, ltl)
%
%   METHODS (opt.Reg)
%     'none'      no action
%     'tikhonov'  Nr = Nll + lambda*L'*L        lambda by GCV if not supplied
%     'tsvd'      truncation (pinv): LSQVCE.ENGINE.BLOCKACCUMULATOR solves on
%                 the eigenvectors of the UNFLOORED Nll above TruncTol*lambda_max
%                 and counts that rank in dof (WP-133) -- of the Jacobi-
%                 equilibrated Nll since WP-151, so the cut does not depend on
%                 the local parameters' units (y is minimum-norm in that
%                 metric; x and v'Pv do not depend on it). The Nr returned here is
%                 Nll floored at TruncTol*lambda_max -- an invertible
%                 placeholder no reader inverts (dg.lambdaUsed is the floor)
%     'lm'        Nr = Nll + lambda*diag(Nll)   Levenberg-Marquardt scaling
%     'auto'      dispatch on dg.status:
%                    'ok'              -> none
%                    'ill-conditioned' -> tikhonov
%                    'rank-deficient'  -> tsvd
%
%   INPUT
%     Nll  nl x nl  local normal matrix
%     nl   nl x 1   local right-hand side
%     dg   struct   from LSQVCE.ENGINE.LOCALDIAGNOSTICS
%     opt  struct   .Reg, .RegLambda, .RegMatrix, .TruncTol
%     mObs 1x1      rows of the system (whitened observations) -- GCV's m
%     ltl  1x1      l'l of the whitened observations -- GCV's residual
%
%   OUTPUT
%     Nr, nr        regularised block and rhs (nr is unchanged by all current
%                   methods; it is returned for symmetry with damped-prior
%                   variants that shift the rhs)
%     dg            with .regUsed and .lambdaUsed filled in
%
%   REGULARISATION IS NOT FREE -- READ THIS
%     Pre-elimination forms Ncc_red = Ncc - Nlc'*inv(Nll)*Nlc. Any bias you
%     introduce into inv(Nll) is transmitted into the COMMON parameters. A
%     regularised local block therefore biases the gravity field, not just the
%     nuisance parameters. This is the opposite of the usual intuition that
%     "the locals get eliminated anyway, so who cares".
%
%     The right fix for a structurally inestimable local parameter is to
%     REMOVE it (or to lengthen the arc, or to add a constraint from an
%     independent source), not to damp it. Use regularisation when the block
%     is merely poorly scaled or weakly determined, not when it is singular by
%     construction.
%
%   CHOOSING THE METHOD
%     tikhonov  A genuine prior: you believe y is small (or L*y is small). Set
%               RegMatrix L to encode that belief -- L = I penalises magnitude,
%               a first-difference L penalises roughness. Use when the block is
%               weakly determined but the parameters are physically meaningful.
%     tsvd      No prior: you simply refuse to invert the null space. Use when
%               the block is rank deficient by construction. The truncated
%               directions are left at zero (the minimum-norm solution) rather
%               than being fitted to noise.
%     lm        Scale-invariant damping. Use when the local parameters have
%               wildly different units and L = I would be meaningless.
%
%   CHOOSING LAMBDA
%     If RegLambda is 0 or empty, Tikhonov calls LSQVCE.ENGINE.GCVLAMBDA (the GCV of
%     the system's local problem, the commons at their a-priori values). GCV is
%     automatic but not infallible -- it can under-regularise when the noise
%     is correlated (which it is, here, after whitening only approximately).
%     For production runs, fix lambda from an L-curve on one representative
%     day and then hold it constant, so the regularisation does not drift
%     between days.
%
%   EXAMPLE
%     opt = struct('Reg','tikhonov', 'RegLambda',0, ...      % 0 -> GCV
%                  'RegMatrix',[], 'TruncTol',1e-10);
%     dg  = lsqvce.engine.localDiagnostics(Nll, opt);
%     [Nr, nr, dg] = lsqvce.engine.applyRegularisation(Nll, nl, dg, opt, mObs, ltl);
%     lsqvce.util.log.msg('%s, lambda = %.3e\n', dg.regUsed, dg.lambdaUsed);
%
%   ACCURACY
%     No reference value is pinned: Nr is Nll + lambda*L'*L or
%     Nll + lambda*diag(Nll) as written, and 'tsvd' floors one symmetric eig
%     (the accumulator then truncates the UNFLOORED Nll instead, WP-133).
%     Pinned by behaviour: every branch gives a finite x on near-collinear
%     local columns (collinearEps 1e-10) and 'auto' leaves a healthy block
%     alone (test_reliability (testAllRegularisationBranchesRun,
%     testAutoRegDispatchesOnStatus)). The GCV lambda is
%     LSQVCE.ENGINE.GCVLAMBDA's (see its ACCURACY). WP-139.
%
%   ERRORS
%     lsqvce:applyReg:badLambda     lambda is not a non-negative scalar
%     lsqvce:applyReg:badRegMatrix  the regularisation matrix has the wrong size
%     lsqvce:applyReg:notPD         the regularised normal matrix is not
%                                   positive definite
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:applyReg:sizeMismatch               Nll is XxX but nl has X elements.RegMatrix
%     lsqvce:applyReg:unknown                    Unknown regularisation 'X'. Valid
%
%   See also LSQVCE.ENGINE.LOCALDIAGNOSTICS, LSQVCE.ENGINE.GCVLAMBDA, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-151 (review 12.84, S1): 'tsvd' is described as the
%   accumulator now cuts it, on the Jacobi-equilibrated Nll; 'auto' reads the
%   equilibrated status of LSQVCE.ENGINE.LOCALDIAGNOSTICS.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: 'tsvd' is described as what
%   it is, the accumulator's truncation of the unfloored Nll (WP-133); the
%   floored Nr returned here is a placeholder no reader inverts. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 finding 16): takes
%   mObs and l'l and hands them to gcvLambda, whose functional is now the
%   GCV of the data (it read N as the design). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(6,6);

    nlp = size(Nll,1);
    if numel(nl) ~= nlp
        error('lsqvce:applyReg:sizeMismatch', ...
              'Nll is %dx%d but nl has %d elements.', nlp, nlp, numel(nl));
    end

    L = lsqvce.util.getfielddef(opt, 'RegMatrix', []);
    if isempty(L)
        L = eye(nlp);
    elseif size(L,2) ~= nlp
        error('lsqvce:applyReg:badRegMatrix', ...
              'RegMatrix has %d columns but the local block has %d parameters.', ...
              size(L,2), nlp);
    end

    tolT = lsqvce.util.getfielddef(opt, 'TruncTol', 1e-10);
    reg  = lower(lsqvce.util.getfielddef(opt, 'Reg', 'none'));

    if strcmp(reg, 'auto')
        switch dg.status
            case 'ok',              reg = 'none';
            case 'ill-conditioned', reg = 'tikhonov';
            otherwise,              reg = 'tsvd';
        end
    end

    Nr = Nll;
    nr = nl;

    switch reg
    case 'none'
        % nothing

    case 'tikhonov'
        lam = lsqvce.util.getfielddef(opt, 'RegLambda', 0);
        if isempty(lam) || all(lam == 0)
            lam = lsqvce.engine.gcvLambda(Nll, nl, L, mObs, ltl);   % WP-135: the data's GCV
        end
        if isscalar(lam)
            validateattributes(lam, {'double'}, {'scalar','finite','nonnegative'}, ...
                               mfilename, 'RegLambda');
            Nr = Nll + lam*(L.'*L);
        else
            if numel(lam) ~= size(L,1)
                error('lsqvce:applyReg:badLambda', ...
                    ['A vector RegLambda must have one entry per row of ' ...
                     'RegMatrix (%d), got %d.'], size(L,1), numel(lam));
            end
            Nr = Nll + L.'*diag(lam(:))*L;
        end
        dg.regUsed    = 'tikhonov';
        dg.lambdaUsed = max(lam);

    case 'lm'
        lam = lsqvce.util.getfielddef(opt, 'RegLambda', 0);
        if isempty(lam) || lam == 0, lam = 1e-8; end
        Nr            = Nll + lam*diag(diag(Nll));
        dg.regUsed    = 'lm';
        dg.lambdaUsed = lam;

    case 'tsvd'
        [V, d] = eig((Nll+Nll.')/2, 'vector');   % symmetric -> V orthonormal
        d      = real(d);
        dm     = max(d);
        if dm <= 0
            error('lsqvce:applyReg:notPD', ...
                  'The local block has no positive eigenvalue (max = %g).', dm);
        end
        fl        = tolT * dm;
        d(d < fl) = fl;                          % floor, keeps Nr invertible
        Nr        = real(V)*diag(d)*real(V).';
        dg.regUsed    = 'tsvd';
        dg.lambdaUsed = fl;

    otherwise
        error('lsqvce:applyReg:unknown', ...
            ['Unknown regularisation ''%s''. Valid: none, tikhonov, tsvd, ' ...
             'lm, auto.'], reg);
    end

    Nr = (Nr + Nr.')/2;
end
