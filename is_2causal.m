function c2_constraints = is_2causal(W, dims, parties, tol)

    if ~exist('tol','var')
        tol = 1e-6;
    end

    if iscell(W)
        W = W{1};
    end

    d = prod(dims);
    actors = parties(2:end-1);
    N = length(parties) - 2;

    c2_constraints = true;

    Wk = cell(1, N);
    S = zeros(d);

    subs = subsets(actors);

    for i = 1:length(subs)
        if isempty(subs{i}) || length(subs{i}) == N
            continue;
        else

            if i == length(subs)
                Wk{1} = W - sum(cat(3, Wk{1:length(subs)-1}), 3);
            else
                Wk{i} = sdpvar(d,d, 'hermitian', 'complex');
            end

            c2_constraints = et(c2_constraints, is_PSD(Wk{i}, tol));
            c2_constraints = et(c2_constraints, is_valid_except_P(Wk{i}, dims, parties, tol));
            c2_constraints = et(c2_constraints, is_compatible(Wk{i}, subs{i}, actorsexclusion(actors, subs{i}), dims, parties, tol));
            S = S + Wk{i};
        end
    end

    c2_constraints = et(c2_constraints, W == S);
    
end