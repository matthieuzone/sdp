function compatible_constraints = is_compatible(W, K1, K2, dims, parties, tol)

    addpath('utils')

    if iscell(W)
        W = W{1};
    end
    
    K1 = wrap_actor(K1);
    K2 = wrap_actor(K2);

    assert(isempty(inter(K1, K2)), 'Error: K1 and K2 must be disjoint sets.');

    actors = parties(2:end-1);
%    P = parties{1};
    F = parties{end};

    compatible_constraints = [];
    subs = subsets(actorsexclusion(actors, K1));
    for i = 1:length(subs)
        if isempty(inter(subs{i},K2))
            continue;
        else
            compatible_constraints = [compatible_constraints, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, [subs{i}, K1]), F], dims)];
        end
    end

    if exist('tol','var')
        compatible_constraints = nullconstraints(compatible_constraints, tol);
    else
        compatible_constraints = nullconstraints(compatible_constraints);
    end


end
