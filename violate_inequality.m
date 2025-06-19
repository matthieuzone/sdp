function [nW, nM, s] = opt(W, M)

    ineq = [5 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0];
    ops = sdpsettings('solver','mosek','verbose',0);

    table = reshape(ineq(2:end), 8, 8);
    score = @(x) sum(x .* table, 'all');

    d = [2, 2, 2, 2, 2, 2];
    parties = {{[]}, {1, 2}, {3, 4}, {5, 6}, {[]}};

    cst = is_2causal(W, d, parties);
    cst = et(cst, nullconstraints(trace(W) - 1));
    for p = 1:3
        for x = 1:2
            cst = et(cst, is_PSD(M{p}{x}{1}));
            cst = et(cst, is_PSD(M{p}{x}{2}));
            cst = et(cst, nullconstraints(PartialTrace(M{p}{x}{1} + M{p}{x}{2}, 1, [2 2]) - eye(2)));
        end
        cst = et(cst, is_PSD(M{p}{2}{1}));
        cst = et(cst, nullconstraints(PartialTrace(M{p}{2}{1}, 1, [2 2]) - eye(2)));
    end
    proba = sdpvar(8,8,'full');
    for o = 0:7
        for i = 0:7
            bo = dec2bin(o, 3) - '0';
            bi = dec2bin(i, 3) - '0';
            proba(i+1,o+1) = real(trace(W * Tensor(M{1}{bi(1)+1}{bo(1)+1}, M{2}{bi(2)+1}{bo(2)+1}, M{3}{bi(3)+1}{bo(3)+1})));
        end
    end
    optimize(cst, score(proba), ops);
    nW = recvalue(W);
    nM = recvalue(M);
    s = score(value(proba));
end

function v = recvalue(x)
    if iscell(x)
        v = cellfun(@recvalue, x, 'UniformOutput', false);
    elseif isa(x, 'sdpvar')
        v = value(x);
    else
        v = x;
    end
end

M = {{rpair(), {rsim(), zeros(4,4)}}, 
     {rpair(), {rsim(), zeros(4,4)}}, 
     {rpair(), {rsim(), zeros(4,4)}}};

for iter = 1:10
    disp(['Iteration: ', num2str(iter)]);

    disp('Optimizing W');
    W = sdpvar(2^6, 2^6, 'hermitian', 'complex');
    [W, M, s] = opt(W, M);
    disp(['Score: ', num2str(s)]);

    for i = 1:3
        for a = 1:2
            disp(['Optimizing M{', num2str(i), '}{1}{', num2str(a), '}']);
            M{i}{1}{a} = sdpvar(4,4);
            [W, M, s] = opt(W, M);
            disp(['Score: ', num2str(s)]);
        end
        disp(['Optimizing M{', num2str(i), '}{2}{1}']);
        M{i}{2}{1} = sdpvar(4,4);
        [W, M, s] = opt(W, M);
        disp(['Score: ', num2str(s)]);
    end
end