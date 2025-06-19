function result = inter(x, y)
    % inter - Intersects indexes of actors in x with those in y.
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

    % Step 2: Intersect each actor in x with y_indexes
    result = {};
    for i = 1:length(x)
        actor_x = x{i};
        intersected_actor = cell(1, length(actor_x));

        total_remaining = 0;
        for j = 1:length(actor_x)
            common = intersect(actor_x{j}, y_indexes);
            intersected_actor{j} = common;
            total_remaining = total_remaining + numel(common);
        end

        if total_remaining > 0
            result{end+1} = intersected_actor; %#ok<AGROW>
        end
    end
end
