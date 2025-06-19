function result = actorsexclusion(x, y)
    % actor_diff_prune - Removes indexes from actors in x if present in y.
    % Handles actors that are 1x1 or 1x2 cell arrays of index lists.
    %
    % Inputs:
    %   x - 1D cell array of actors
    %   y - 1D cell array of actors
    %
    % Output:
    %   result - pruned cell array of actors

    % Step 1: Gather all indexes from y

    x = wrap_actor(x);
    y = wrap_actor(y);

    y_indexes = [];
    for i = 1:length(y)
        actor_y = y{i};
        for j = 1:length(actor_y)
            y_indexes = [y_indexes, actor_y{j}];
        end
    end
    y_indexes = unique(y_indexes);

    % Step 2: Prune each actor in x
    result = {};
    for i = 1:length(x)
        actor_x = x{i};
        pruned_actor = cell(1, length(actor_x));

        total_remaining = 0;
        for j = 1:length(actor_x)
            pruned_indexes = setdiff(actor_x{j}, y_indexes);
            pruned_actor{j} = pruned_indexes;
            total_remaining = total_remaining + numel(pruned_indexes);
        end

        if total_remaining > 0
            result{end+1} = pruned_actor; %#ok<AGROW>
        end
    end
end


