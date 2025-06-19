function rpair = rpair()
    good = false;
    while ~good
        A = rpsd(4);
        B = rpsd(4);
        B = B - tr_replace(A + B, 1, [2 2]) + eye(4)/2;
        rpair = {A, B};
        good = is_PSD(rpair{1}) && is_PSD(rpair{2});
    end
end

function A = rpsd(n)
    A = randn(n);
    A = A*A';
    A = A/(2*trace(A));
end
