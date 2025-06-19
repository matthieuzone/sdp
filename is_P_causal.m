function P_constraints = is_P_causal(W, P, dims, parties, tol)

    if ~exist('tol','var')
        tol = 1e-6;
    end
    if iscell(W)
        W = W{1};
    end

    
    actors = parties(2:end-1);
    F = parties{end};
    
    d = prod(dims);
    N = length(P);
    P_constraints = true;

    Wk = cell(1, N);

    for i = 1:N
        if i == N
            Wk{1} = W - sum(cat(3, Wk{1:N-1}), 3);
        else
            Wk{i} = sdpvar(d,d, 'hermitian', 'complex');
        end

        P_constraints = et(P_constraints, is_PSD(Wk{i}, tol));
        P_constraints = et(P_constraints, is_valid(Wk{i}, dims, parties),tol);
        P_constraints = et(P_constraints, is_compatible(Wk{i}, P{i}, [actorsexclusion(actors, P{i}), F], dims, parties, tol));
        %pour toute matrice, machin est causal
        %pour tout rho
    end

    assign(Wk{1}, W);
    assign(Wk{2}, 0);

end
