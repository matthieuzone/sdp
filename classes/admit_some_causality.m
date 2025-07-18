function sc_constraints = admit_some_causality(W, dims, parties, tol)
    %admit_some_causality - Create constraints for admitting some causality in a superoperator W
    %even if W is a numerical value, it will return the constraints for the convex combination without solving them.

    if ~exist('tol','var')
        tol = 1e-6;
    end

    if iscell(W)
        W = W{1};
    end

    d = prod(dims);
    actors = parties(2:end-1); %exclude P and F
    N = length(parties) - 2;

    sc_constraints = true;

    S = zeros(d);

    subs = subsets(actors);
    Wk = cell(1, length(subs));

    %For each K_1 and K_2, add a term compatible with K_1 < K_2 (not valid but with the condition except for P) to the convex combination
    for i = 1:length(subs)

        if isempty(subs{i}) || length(subs{i}) == N
            continue;
        else

            NmK =  actorsexclusion(actors, subs{i});
            K2s = subsets(NmK);
            Wk{i} = cell(1, length(K2s));

            for j = 1:length(K2s) 
                if isempty(K2s{j})
                    continue;
                end
                Wk{i}{j} = sdpvar(d,d, 'hermitian', 'complex');
                sc_constraints = et(sc_constraints, is_valid_except_P(Wk{i}{j}, dims, parties, tol));
                sc_constraints = et(sc_constraints, is_compatible(Wk{i}{j}, subs{i}, K2s{j}, dims, parties, tol));
                S = S + Wk{i}{j}; 
            end
        end
    end

    sc_constraints = et(sc_constraints, W == S);
    
end