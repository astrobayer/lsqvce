function [Aw, Bw, lw] = activeColumns(whiten, A, B, l)
%ACTIVECOLUMNS  Run a dense whitener on the non-zero columns of sparse blocks.
%   [Aw,Bw,lw] = activeColumns(whiten, A, B, l)
%
%   The ONE sparse path of the kernel whiteners (WP-123): LSQVCE.KERNEL.
%   SEPARABLEWHITEN, LSQVCE.KERNEL.VARWHITEN and LSQVCE.KERNEL.ARMAWHITEN
%   reshape to N-D arrays or call filter, which MATLAB refuses for sparse
%   input (MATLAB:reshape:sparseNDmatrix) -- yet LSQVCE.ENGINE.MAKELOADER
%   promises that B may be sparse, and LSQVCE.ENGINE.SYNTHLOADER recommends
%   it. Each whitener is a ROW operator applied column by column, so a
%   column that is zero in the batch stays zero: the dense whitener runs on
%   the columns with a non-zero entry only, and the result is scattered back
%   into a sparse matrix of the original width. Memory is batch rows x
%   active columns, not batch rows x all locals of the system.
%
%   INPUT
%     whiten   1x1 function_handle  @(A,B,l) -> [Aw,Bw,lw], the whitener on
%              DENSE blocks with its model arguments bound
%     A        m x nc  double  common design rows (dense or sparse)
%     B        m x nl  double  local design rows  (dense or sparse, may be [])
%     l        m x 1   double  observations       (dense or sparse)
%
%   OUTPUT
%     Aw, Bw, lw   the whitened blocks; each is sparse exactly when its
%                  input was, with the input's column count
%
%   ERRORS  none of its own (the whitener's errors pass through).
%
%   EXAMPLE
%     W = @(A,B,l) deal(2*A, 2*B, 2*l);      % stand-in for a whitener
%     [~, Bw] = activeColumns(W, zeros(4,0), sparse([1 0; 0 0; 0 0; 0 0]), ones(4,1));
%     full(Bw)                               % [2 0; 0 0; 0 0; 0 0], sparse
%
%   See also LSQVCE.KERNEL.SEPARABLEWHITEN, LSQVCE.KERNEL.VARWHITEN,
%            LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.ENGINE.MAKELOADER
%
%   Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-123.

    [Ac, ca] = compress(A);
    [Bc, cb] = compress(B);
    [Awc, Bwc, lwc] = whiten(Ac, Bc, full(l));
    Aw = expand(Awc, ca, A);
    Bw = expand(Bwc, cb, B);
    lw = lwc;
    if issparse(l), lw = sparse(lwc); end
end

function [Xc, c] = compress(X)
%COMPRESS  The dense non-zero columns of a sparse X (X itself if dense).
    c = [];
    Xc = X;
    if issparse(X)
        c = find(any(X, 1));
        Xc = full(X(:, c));
    end
end

function Xw = expand(Xwc, c, X)
%EXPAND  Scatter whitened active columns back to X's sparse layout.
    Xw = Xwc;
    if issparse(X)
        [i, j, v] = find(Xwc);
        cj = c(j);
        Xw = sparse(i(:), cj(:), v(:), size(Xwc, 1), size(X, 2));
    end
end
