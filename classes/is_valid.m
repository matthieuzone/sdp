function valid_constraints = is_valid(W, dims, parties, tol)

    if iscell(W)
        W = W{1};
    end

    valid_constraints = [];

    actors = parties(2:end-1);
    P = {[],parties{1}{1}};
    F = parties{end};

    subs = subsets(actors);
    for i = 1:length(subs)
        if isempty(subs{i})
            continue;
        else
            valid_constraints = [valid_constraints, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, subs{i}), F], dims)];
        end
    end

    valid_constraints = [valid_constraints, reduced_for_constraint(W, P, [actors, F], dims)];

    if exist('tol','var')
        valid_constraints = nullconstraints(valid_constraints, tol);
        valid_constraints = et(valid_constraints, is_PSD(W, tol));
    else
        valid_constraints = nullconstraints(valid_constraints);
        valid_constraints = et(valid_constraints, is_PSD(W));
    end

end