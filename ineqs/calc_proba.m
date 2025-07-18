function proba = calc_proba(W, M)
    %calc_proba - Calculate the probability distribution from a superoperator W and instruments M.

    proba = sdpvar(8,8,'full');
    for o = 0:7
        for i = 0:7
            bo = dec2bin(o, 3) - '0';
            bi = dec2bin(i, 3) - '0';
            proba(i+1,o+1) = real(trace(W * Tensor(M{1}{bi(1)+1}{bo(1)+1}, M{2}{bi(2)+1}{bo(2)+1}, M{3}{bi(3)+1}{bo(3)+1})));
        end
    end