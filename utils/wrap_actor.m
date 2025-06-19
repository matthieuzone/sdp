function out = wrap_actor(x)
    % wrap_actor - Ensures x is a cell array of actors
    % An actor is a 1x1 or 1x2 cell array of index lists
    %
    % Input:
    %   x - either a single actor or a cell array of actors
    % Output:
    %   out - always a cell array of actors

    if iscell(x) && all(cellfun(@iscell, x))
        % Already a cell array of actors
        out = x;
    else
        % Wrap single actor in a cell array
        out = {x};
    end
end