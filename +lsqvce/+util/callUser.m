function [nOut, varargout] = callUser(fh, minOut, maxOut, ctx, varargin)
%CALLUSER  Call a user-supplied function: outputs counted ONCE, errors wrapped.
%   [nOut, out1, ..., outMax] = lsqvce.util.callUser(fh, minOut, maxOut, ctx, args...)
%
%   THE one place the toolbox calls a function the user wrote (batch sources,
%   plot hooks, arc resolvers; buildDay goes through makeLoader, which applies
%   the same two rules). The two rules, R5 error doctrine (theory_r5 ch. 7,
%   finding H-08):
%     1. HOW MANY OUTPUTS the function returns is asked, not probed on every
%        call. nargout(fh) answers for named functions without varargout; it
%        answers -1 for anonymous handles and varargout functions, and THEN --
%        and only then -- the count is found by ONE trial call that asks for
%        maxOut outputs and steps down on exactly the three MATLAB "too many
%        outputs" identifiers. The call's outputs are RETURNED, so the trial is
%        never wasted, and the caller stores nOut and passes minOut == maxOut
%        == nOut from then on (no probing ever again).
%     2. ANY OTHER exception inside the user's function is re-raised as
%        lsqvce:callUser:threw with the call context in the message and the
%        user's exception attached UNMODIFIED as the cause (ME.cause{1}), so
%        the user sees both where the toolbox was and what their code said.
%
%   INPUT
%     fh      1x1 function_handle   the user's function
%     minOut  1x1 double >= 0       fewest outputs the caller can work with
%     maxOut  1x1 double >= minOut  most outputs the caller wants
%     ctx     char/string           what is being called, for messages
%                                   (e.g. 'batchFcn(k)', 'plot hook after/kbr')
%     args    any                   the arguments of the call
%   OUTPUT
%     nOut         1x1 double   how many outputs fh actually returned
%                               (minOut <= nOut <= maxOut)
%     out1..outMax any          the outputs of THIS call; positions above nOut
%                               are [] (the caller asked for maxOut names)
%
%   ERRORS
%     lsqvce:callUser:badHandle   fh is not a function handle
%     lsqvce:callUser:tooFew      fh cannot return minOut outputs (nargout(fh)
%                                 is known and smaller, or every trial down to
%                                 minOut raised a too-many-outputs identifier)
%     lsqvce:callUser:threw       fh raised anything else; the original is the
%                                 cause
%
%   ACCURACY  none of its own: outputs pass through untouched.
%
%   EXAMPLE
%     % first batch: learn n; later batches: fixed n, no probing
%     [n, A, l, p] = lsqvce.util.callUser(batchFcn, 2, 3, 'batchFcn(k)', 1);
%     for k = 2:nBatch
%         [~, A, l, p] = lsqvce.util.callUser(batchFcn, n, n, 'batchFcn(k)', k);
%     end
%
%   See also NARGOUT, LSQVCE.ENGINE.MAKELOADER
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M1, WP-80.

    if ~isa(fh, 'function_handle')
        error('lsqvce:callUser:badHandle', ...
              '%s: expected a function handle, got %s.', char(ctx), class(fh));
    end
    minOut = double(minOut);  maxOut = double(maxOut);
    varargout = cell(1, maxOut);

    % Rule 1a: ask.
    declared = nargout(fh);                     % -1: anonymous / varargout
    if declared >= 0
        if declared < minOut
            error('lsqvce:callUser:tooFew', ...
                  '%s: the function returns %d output(s); at least %d are needed.', ...
                  char(ctx), declared, minOut);
        end
        n = min(declared, maxOut);
        [varargout{1:n}] = guarded(fh, ctx, varargin{:});
        nOut = n;
        return
    end

    % Rule 1b: probe ONCE, stepping down on the three arity identifiers only.
    arity = {'MATLAB:TooManyOutputs', 'MATLAB:maxlhs', 'MATLAB:unassignedOutputs'};
    for n = maxOut:-1:minOut
        try
            if n == 0
                fh(varargin{:});
            else
                [varargout{1:n}] = fh(varargin{:});
            end
            nOut = n;
            return
        catch ME
            if n > minOut && any(strcmp(ME.identifier, arity))
                continue                        % fewer outputs, once
            elseif any(strcmp(ME.identifier, arity))
                error('lsqvce:callUser:tooFew', ...
                      '%s: the function cannot return %d output(s) (%s).', ...
                      char(ctx), minOut, ME.identifier);
            else
                rethrowWrapped(ME, ctx);
            end
        end
    end
end

function varargout = guarded(fh, ctx, varargin)
%GUARDED  Call fh with the caller's nargout; wrap anything it throws.
    try
        if nargout == 0
            fh(varargin{:});
        else
            [varargout{1:nargout}] = fh(varargin{:});
        end
    catch ME
        rethrowWrapped(ME, ctx);
    end
end

function rethrowWrapped(ME, ctx)
%RETHROWWRAPPED  Rule 2: lsqvce:callUser:threw with the user's error as cause.
    E = MException('lsqvce:callUser:threw', ...
                   '%s threw %s: %s', char(ctx), ME.identifier, ME.message);
    E = addCause(E, ME);
    throw(E);
end
