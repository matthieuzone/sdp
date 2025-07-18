function compatible_constraints = is_compatible(W, K1, K2, dims, parties, tol)

    %is_compatible - Create constraints for a superoperator W to be compatible with K1 and K2.
    %If W is a numerical value, it will return a boolean

    if iscell(W)
        W = W{1};
    end
    
    %ensure K1 and K2 are cell arrays of actors
    K1 = wrap_actor(K1);
    K2 = wrap_actor(K2);

    assert(isempty(inter(K1, K2)), 'Error: K1 and K2 must be disjoint sets.');

    actors = parties(2:end-1);
    F = parties{end};

    ctr = [];
    subs = subsets(actorsexclusion(actors, K1));
    for i = 1:length(subs) %forall X in N such that X intersect K_2
        if isempty(inter(subs{i},K2))
            continue;
        else
            ctr = [ctr, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, [subs{i}, K1]), F], dims)]; %[1 - AO]^X [AIO]^(N\K_1\X U F) W
        end
    end

    if exist('tol','var')
        compatible_constraints = nullconstraints(ctr, tol); % = 0
    else
        compatible_constraints = nullconstraints(ctr);
    end


end
