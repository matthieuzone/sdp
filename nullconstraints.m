function constraints = nullconstraints(constraints, tol)

    if ~exist('tol','var')
        tol = 1e-6;
    end

    if isa(constraints, 'sdpvar')
        constraints = constraints == 0; 
    else
        constraints = max(abs(constraints),[],"all") <= tol;
    end
end