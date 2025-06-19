function score = score(W, M, ineq)

    table = reshape(ineq(2:end), 8, 8);
    score = sum(proba(W,M) .* table, 'all');
end