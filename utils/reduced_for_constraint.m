function reduced = reduced_for_constraint(W, AOs, AIOs, dims)
    % This function returns [1 - AOs][AIOs]W
    % AOs: cell array of actors who's outputs are to be traced out
    % AIOs: cell array of actors who are to be traced out

    AOs = wrap_actor(AOs);
    AIOs = wrap_actor(AIOs);

    AIOs = cell2list(AIOs);
    W = tr_replace(W, AIOs, dims);
    reduced = removeAO(W, AOs, dims);
end

function removedAO = removeAO(W, X, dims)
    removedAO = W;
    for k = 1:length(X)
        removedAO = removedAO - tr_replace(removedAO, X{k}{2}, dims);
    end
end

function list = cell2list(cellArray)
    % Convert a (possibly nested) cell array to a list of elements
    list = [];
    for i = 1:numel(cellArray)
        if iscell(cellArray{i})
            list = [list, cell2list(cellArray{i})]; % Recursively handle nested cells
        else
            list = [list, cellArray{i}];
        end
    end
end

function prettyPrintCellArray(cellArray, indent)
    if nargin < 2
        indent = 0; % Default indentation level
    end
    for i = 1:numel(cellArray)
        if iscell(cellArray{i})
            fprintf('%sCell {\n', repmat(' ', 1, indent));
            prettyPrintCellArray(cellArray{i}, indent + 4); % Increase indentation for nested cells
            fprintf('%s}\n', repmat(' ', 1, indent));
        else
            fprintf('%s%s\n', repmat(' ', 1, indent), mat2str(cellArray{i}));
        end
    end
end