ineqs = readmatrix('data/3lazy_2-causal_repr.txt');

for i = 1:size(ineqs, 1)
    ineq = ineqs(i, :);
    disp(['Inequality ', num2str(i), ':']);
    [W, M] = max_inequality(ineq, '2-causal', 1);
    s = score(calc_proba(W, M), ineq);
    disp(['Violated: ']);
    disp(s < -ineq(1) - 1e-3);
end