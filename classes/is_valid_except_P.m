function valid_constraints = is_valid_except_P(W, dims, parties, tol)
    
    valid_constraints = [];

    actors = parties(2:end-1);
    F = parties{end};
    if iscell(W)
        W = W{1};
    end

    subs = subsets(actors);
    for i = 1:length(subs)
        if isempty(subs{i})
            continue;
        else
            valid_constraints = [valid_constraints, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, subs{i}), F], dims)];
        end
    end

    if exist('tol','var')
        valid_constraints = nullconstraints(valid_constraints, tol);
        valid_constraints = et(valid_constraints, is_PSD(W, tol));
    else
        valid_constraints = nullconstraints(valid_constraints);
        valid_constraints = et(valid_constraints, is_PSD(W));
    end
end