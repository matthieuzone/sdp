function constraints = nullconstraints(eqs, tol)
    %create the constraints eqs == 0 for sdpvars or retrun the truth value for numerical values

    if ~exist('tol','var')
        tol = 1e-6;
    end

    if isa(eqs, 'sdpvar')
        constraints = eqs == 0; 
    else
        constraints = max(abs(eqs),[],"all") <= tol;
    end
end