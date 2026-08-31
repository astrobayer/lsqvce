function blk = leverageBlocks(Ncc_red, Nll, Nlc)
%LSQVCE.ENGINE.LEVERAGEBLOCKS  Precompute the small matrices LSQVCE.ENGINE.LEVERAGE needs.
%
%   blk = LSQVCE.ENGINE.LEVERAGEBLOCKS(Ncc_red, Nll, Nlc)
%
%   Called ONCE PER SYSTEM, then reused for every arc of that system. All of the
%   O(nc^3) and O(nl^3) work happens here; LSQVCE.ENGINE.LEVERAGE is then O(m*nc^2/2) per
%   arc with no further factorisations.
%
%   INPUT
%     Ncc_red  nc x nc  the TOTAL reduced normal matrix, summed over ALL systems
%                       (this is the joint-inverse part -- do not pass one
%                       system's Ncc)
%     Nll      nl x nl  this system's local normal matrix (regularised, as used)
%     Nlc      nl x nc  this system's cross block
%
%   OUTPUT
%     blk.Lc   chol(Ncc_red, 'lower')
%     blk.Ll   chol(Nll,     'lower')
%     blk.G    inv(Nll) * Nlc            nl x nc
%     blk.SG   Sxx * G'                  nc x nl
%     blk.GSG  G * Sxx * G'              nl x nl
%
%   STALENESS
%     These blocks depend on the variance components through the whitening. They
%     MUST be rebuilt at every VCE iteration. LSQVCE.ENGINE.LEVERAGE warns if it detects a
%     mismatch (leverage outside [0,1]), but it cannot always tell.
%
%   EXAMPLE
%     blk = lsqvce.engine.leverageBlocks(Ncc_red_total, Nll_s, Nlc_s);
%     for each arc of system s
%         h = lsqvce.engine.leverage(Aw, Bw, blk);
%     end
%
%   ERRORS
%     lsqvce:leverageBlocks:sizeMismatch  the blocks do not line up
%     lsqvce:leverageBlocks:cNotPD        the common block is not positive definite
%     lsqvce:leverageBlocks:lNotPD        a local block is not positive definite
%
%   See also LSQVCE.ENGINE.LEVERAGE, LSQVCE.ENGINE.REDUNDANCY
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(3,3);

    nc = size(Ncc_red,1);
    nl = size(Nll,1);
    if ~isequal(size(Nlc), [nl nc])
        error('lsqvce:leverageBlocks:sizeMismatch', ...
              'Nlc must be %d x %d, got %d x %d.', nl, nc, size(Nlc,1), size(Nlc,2));
    end

    [Lc, p1] = chol((Ncc_red + Ncc_red.')/2, 'lower');
    if p1 > 0
        error('lsqvce:leverageBlocks:cNotPD', ...
            ['The reduced common normal matrix is not positive definite. The ' ...
             'common block has a datum defect -- add constraints, or supply a ' ...
             'prior via the common-part regularisation.']);
    end

    [Ll, p2] = chol((Nll + Nll.')/2, 'lower');
    if p2 > 0
        error('lsqvce:leverageBlocks:lNotPD', ...
            ['The local normal matrix is not positive definite. Pass the ' ...
             'REGULARISED Nll (the one actually used in pre-elimination), not ' ...
             'the raw one.']);
    end

    G   = Nll \ Nlc;                       % nl x nc
    SG  = Lc.' \ (Lc \ G.');               % nc x nl  = Sxx * G'
    GSG = G * SG;                          % nl x nl  = G * Sxx * G'
    GSG = (GSG + GSG.')/2;

    blk = struct('Lc', Lc, 'Ll', Ll, 'G', G, 'SG', SG, 'GSG', GSG);
end
