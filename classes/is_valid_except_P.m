function valid_constraints = is_valid_except_P(W, dims, parties, tol)
    %is_valid_except_P - Create constraints for a superoperator W to be valid except for P, used in the characterization of 2-causal and subsets-causal processes. returns a boolean if W is a numerical value.

    ctr = [];

    actors = parties(2:end-1);
    F = parties{end};
    if iscell(W)
        W = W{1};
    end

    subs = subsets(actors);
    for i = 1:length(subs) % For all X in N
        if isempty(subs{i})
            continue;
        else
            % [1 - AO]^X [AIO]^(N\X U F) W
            ctr = [ctr, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, subs{i}), F], dims)];
        end
    end

    % = 0
    if exist('tol','var')
        valid_constraints = nullconstraints(ctr, tol);
        valid_constraints = et(valid_constraints, is_PSD(W, tol));
    else
        valid_constraints = nullconstraints(ctr);
        valid_constraints = et(valid_constraints, is_PSD(W));
    end
end