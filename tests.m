ineqs = readmatrix('data/3lazy_2-causal_repr.txt');
t = 0;

for i = 1:size(ineqs, 1)
    ineq = ineqs(i, :);
    disp(['Inequality ', num2str(i), ':']);
    [W, M] = max_inequality(ineq, 'valid', 2, true);
    s = score(calc_proba(W, M), ineq);
    disp('Violated: ');
    disp(s < -ineq(1) - 1e-3);
    if s < -ineq(1) - 1e-3
        t = t + 1;
    end
end
disp(['Total violated inequalities: ', num2str(t)]);