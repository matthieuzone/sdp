function valid_constraints = is_valid(W, dims, parties, tol)
    %is_valid - Create constraints for a superoperator W to be valid, returns a boolean if W is a numerical value.
    if iscell(W)
        W = W{1};
    end

    ctr = [];

    actors = parties(2:end-1);
    P = {[],parties{1}{1}};
    F = parties{end};

    subs = subsets(actors); 
    for i = 1:length(subs) % For all X in N
        if isempty(subs{i})
            continue;
        else
            % [1 - AO]^X [AIO]^(N\X U F) W
            ctr = [ctr, reduced_for_constraint(W, subs{i}, [actorsexclusion(actors, subs{i}), F], dims)];
        end
    end

    % [1 - AO]^P [AIO]^(N U F) W
    ctr = [ctr, reduced_for_constraint(W, P, [actors, F], dims)];

    % = 0
    if exist('tol','var')
        valid_constraints = nullconstraints(ctr, tol);
        valid_constraints = et(valid_constraints, is_PSD(W, tol));
    else
        valid_constraints = nullconstraints(ctr);
        valid_constraints = et(valid_constraints, is_PSD(W));
    end

end