function [W, M] = max_inequality(ineq, class, niter)

    ops = sdpsettings('solver','mosek','verbose',0);

    d = [2, 2, 2, 2, 2, 2];
    p = {{[]}, {1, 2}, {3, 4}, {5, 6}, {[]}};

    
    idd = [1,0,0,1];
    id_CJ = pure_to_mixed(idd);
    p0 = pure_to_mixed([1,0]);
    p1 = pure_to_mixed([0,1]);
    t01 = Tensor(p0, p1);
    M = {{{t01, t01}, {id_CJ, zeros(4,4)}}, 
        {{t01, t01} , {id_CJ, zeros(4,4)}}, 
        {{t01, t01}, {id_CJ, zeros(4,4)}}};
    % M = {{{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
    %     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
    %     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}};

    for iter = 1:niter
        disp(['Iteration: ', num2str(iter)]);

        %disp('Optimizing W');
        W = sdpvar(2^6, 2^6, 'hermitian', 'complex');
        switch class
            case 'causal'
                cst = superop_in_QCCC_cone(W, d, p);
            case 'QCCC'
                cst = superop_in_QCCC_cone(W, d, p);
            case '2-causal'
                cst = is_2causal(W, d, p);
            case 'some-causality'
                cst = admit_some_causality(W, d, p);
            case 'valid'
                %cst = superop_in_valid_cone(W, d, p);
                cst = is_valid(W, d, p);
            otherwise
                error('Unknown class');
        end
        cst = et(cst, nullconstraints(trace(W) - 8));

        optimize(cst, score(calc_proba(W, M), ineq), ops);
        W = value(W);
        %disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);

        for i = 1:3
            %disp(['Optimizing M{', num2str(i), '}{1}']);
            M{i}{1}{1} = sdpvar(4,4, 'hermitian', 'complex');
            M{i}{1}{2} = sdpvar(4,4, 'hermitian', 'complex');
            cst = et(is_PSD(M{i}{1}{2}), is_PSD(M{i}{1}{1}));
            cst = et(cst, nullconstraints(PartialTrace(M{i}{1}{1} + M{i}{1}{2}, 2, [2 2]) - eye(2)));
            optimize(cst, score(calc_proba(W, M), ineq), ops);
            M{i}{1}{1} = value(M{i}{1}{1});
            M{i}{1}{2} = value(M{i}{1}{2});
            %disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);

            %disp(['Optimizing M{', num2str(i), '}{2}']);
            M{i}{2}{1} = sdpvar(4,4, 'hermitian', 'complex');
            cst = is_PSD(M{i}{2}{1});
            cst = et(cst, nullconstraints(PartialTrace(M{i}{2}{1}, 2, [2 2]) - eye(2)));
            optimize(cst, score(calc_proba(W, M), ineq), ops);
            M{i}{2}{1} = value(M{i}{2}{1});
        end
        disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
    end
end