function A = rsim()
    good = false;
    while ~good
        A = rpsd(4);
        A = A - tr_replace(A, 1, [2 2]) + eye(4)/2;
        good = is_PSD(A);
    end
end

function A = rpsd(n)
    A = randn(n) + 1i*randn(n);
    A = A*A';
    A = A/(2*trace(A));
end