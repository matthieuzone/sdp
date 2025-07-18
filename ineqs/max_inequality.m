function [W, M] = max_inequality(ineq, class, niter, optimize_instruments, instruments_init, d, verbose)
    %max_inequality - Optimize a superoperator W and instruments M to maximize the score against an inequality.
    %    % Inputs:
    %    ineq - The inequality we try to violate, given as a vector with first the bound and then the coefficients.
    %    class - The class of processes we use to try and violate it ('causal', 'QCCC', '2-causal', 'some-causality', 'valid').
    %    niter - Number of iterations for the optimization.
    %    optimize_instruments - Whether to optimize the instruments M or only W, if false, niter > 1 is useless.
    %    instruments_init - Initial instruments to use, if empty, default ones are used.
    %    d - Dimensions of the systems, if empty, default is [2, 2, 2, 2, 2, 2].
    %    verbose - Whether to display the optimization process.
    ops = sdpsettings('solver','scs','verbose',1);

    if nargin < 7
        verbose = false;
    end
    if nargin < 6
        d = [2, 2, 2, 2, 2, 2];
    end
    if nargin < 5
        % best found instruments in dim 2*2
        idd = [1,0,0,1];
        id_CJ = pure_to_mixed(idd);
        p0 = pure_to_mixed([1,0]);
        p1 = pure_to_mixed([0,1]);
        t01 = Tensor(p0, p1);
        t10 = Tensor(p1, p0);
        instruments_init = {{{id_CJ, zeros(4,4)}, {t01, t10}}, 
            {{id_CJ, zeros(4,4)}, {t01, t10}}, 
            {{id_CJ, zeros(4,4)}, {t01, t10}}};
    end

    p = {{[]}, {1, 2}, {3, 4}, {5, 6}, {[]}};

    M = instruments_init;

    for iter = 1:niter
        if verbose
            disp(['Iteration: ', num2str(iter)]);
        end

        if verbose
            disp('Optimizing W');
        end
        W = sdpvar(prod(d), prod(d), 'hermitian'); %, 'hermitian', 'complex');
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
                cst = superop_in_valid_cone(W, d, p);
                %cst = is_valid(W, d, p);
            otherwise
                error('Unknown class');
        end
        cst = et(cst, nullconstraints(trace(W) - 8));

        optimize(cst, score(calc_proba(W, M), ineq), ops);
        W = value(W);
        if verbose
            disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
        end

        if optimize_instruments
            for i = 1:3
                if verbose
                    disp(['Optimizing M{', num2str(i), '}{1}{1}']);
                end
                M{i}{1}{1} = sdpvar(4,4, 'hermitian', 'complex');
                cst = is_PSD(M{i}{1}{1});
                cst = et(cst, nullconstraints(PartialTrace(M{i}{1}{1}, 2, [2 2]) - eye(2)));
                optimize(cst, score(calc_proba(W, M), ineq), ops);
                M{i}{1}{1} = value(M{i}{1}{1});
                if verbose
                    disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
                end

                if verbose
                    disp(['Optimizing M{', num2str(i), '}{2}{1}']);
                end
                M{i}{2}{1} = sdpvar(4,4, 'hermitian', 'complex');
                cst = is_PSD(M{i}{2}{1});
                cst = et(cst, nullconstraints(PartialTrace(M{i}{2}{1} + M{i}{2}{2}, 2, [2 2]) - eye(2)));
                optimize(cst, score(calc_proba(W, M), ineq), ops);
                M{i}{2}{1} = value(M{i}{2}{1});
                if verbose
                    disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
                end

                if verbose
                    disp(['Optimizing M{', num2str(i), '}{2}{2}']);
                end
                M{i}{2}{2} = sdpvar(4,4, 'hermitian', 'complex');
                cst = is_PSD(M{i}{2}{2});
                cst = et(cst, nullconstraints(PartialTrace(M{i}{2}{1} + M{i}{2}{2}, 2, [2 2]) - eye(2)));
                optimize(cst, score(calc_proba(W, M), ineq), ops);
                M{i}{2}{2} = value(M{i}{2}{2});
                if verbose
                    disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
                end
            end
        end
        if verbose
            disp(['Score: ', num2str(score(calc_proba(W, M), ineq))]);
        end
    end
end