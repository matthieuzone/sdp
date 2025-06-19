function causal_constraints = is_causal(W, dims, p, tol)

    if ~exist('tol','var')
        tol = 1e-6;
    end
    if iscell(W)
        W = W{1};
    end

    addpath('utils')
    d = prod(dims);
    actors = p(2:end-1);
    N = length(p) - 2;

    causal_constraints = true;
     
    Wk = cell(1, N);
    S = zeros(d);

    for i = 1:N
        if i == N
            Wk{1} = W - sum(cat(3, Wk{1:N-1}), 3);
        else
            Wk{i} = sdpvar(d,d, 'hermitian', 'complex');
        end

        causal_constraints = et(causal_constraints, is_PSD(Wk{i}, tol));
        causal_constraints = et(causal_constraints, is_valid(Wk{i}, dims, p, tol));
        causal_constraints = et(causal_constraints, is_compatible(Wk{i}, p{i}, actorsexclusion(actors, p{i}), dims, p, tol));
        S = S + Wk{i};
        %pour toute matrice, machin est causal
        %pour tout rho
    end
    causal_constraints = et(causal_constraints, W == S);

%    assign(Wk{1}, W);
%    assign(Wk{2}, 0);

end