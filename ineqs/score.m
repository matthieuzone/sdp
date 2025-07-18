function score = score(proba, ineq)
    % Calculate the score of a probability distribution against an inequality.

    table = transpose(reshape(ineq(2:end), 8, 8));
    score = sum(proba .* table, 'all');
end