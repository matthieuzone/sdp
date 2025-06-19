function max_inequality(ineq, class, niter)
    ineq = [5 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0];

    ops = sdpsettings('solver','mosek','verbose',0);

    d = [2, 2, 2, 2, 2, 2];
    parties = {{[]}, {1, 2}, {3, 4}, {5, 6}, {[]}};

    table = reshape(ineq(2:end), 8, 8);

    M = {{rpair(), {eye(4)/2, zeros(4,4)}}, 
        {rpair(), {eye(4)/2, zeros(4,4)}}, 
        {rpair(), {eye(4)/2, zeros(4,4)}}};

    for iter = 1:niter
        disp(['Iteration: ', num2str(iter)]);

        disp('Optimizing W');
        W = sdpvar(2^6, 2^6, 'hermitian', 'complex');
        cst = is_2causal(W, d, parties);
        cst = et(cst, nullconstraints(trace(W) - 1));

        optimize(cst, score(W, M, table), ops);
        W = value(W);
        disp(['Score: ', num2str(score(W, M, table))]);

        for i = 1:3
            for a = 1:2
                disp(['Optimizing M{', num2str(i), '}{1}{', num2str(a), '}']);
                M{i}{1}{a} = sdpvar(4,4, 'hermitian', 'complex');
                cst = is_PSD(M{i}{1}{a});
                cst = et(cst, nullconstraints(PartialTrace(M{i}{1}{1} + M{i}{1}{2}, 1, [2 2]) - eye(2)));
                
                optimize(cst, score(W, M, table), ops);
                M{i}{1}{a} = value(M{i}{1}{a});
                disp(['Score: ', num2str(score(W, M, table))]);
            end
        end
    end
end

function score = score(W, M, table)
    proba = sdpvar(8,8,'full');
    for o = 0:7
        for i = 0:7
            bo = dec2bin(o, 3) - '0';
            bi = dec2bin(i, 3) - '0';
            proba(i+1,o+1) = real(trace(W * Tensor(M{1}{bi(1)+1}{bo(1)+1}, M{2}{bi(2)+1}{bo(2)+1}, M{3}{bi(3)+1}{bo(3)+1})));
        end
    end
    score = sum(proba .* table, 'all');
end