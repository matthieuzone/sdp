function score = score(proba, ineq)

    table = transpose(reshape(ineq(2:end), 8, 8));
    score = sum(proba .* table, 'all');
end