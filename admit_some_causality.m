function sc_constraints = admit_some_causality(W, dims, parties, tol)

    if ~exist('tol','var')
        tol = 1e-6;
    end

    if iscell(W)
        W = W{1};
    end

    addpath('utils')
    d = prod(dims);
    actors = parties(2:end-1);
    N = length(parties) - 2;

    sc_constraints = true;

    S = zeros(d);

    subs = subsets(actors);
    Wk = cell(1, length(subs));

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
                sc_constraints = et(sc_constraints, is_PSD(Wk{i}{j}, tol));
                sc_constraints = et(sc_constraints, is_valid_except_P(Wk{i}{j}, dims, parties, tol));
                sc_constraints = et(sc_constraints, is_compatible(Wk{i}{j}, subs{i}, K2s{j}, dims, parties, tol));
                S = S + Wk{i}{j};
            end
        end
    end

    sc_constraints = et(sc_constraints, W == S);
    
end