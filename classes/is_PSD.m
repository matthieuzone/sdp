function PSD_constraints = is_PSD(W, tol)
    
    if ~exist('tol','var')
        tol = 1e-6;
    end
    if iscell(W)
        W = W{1};
    end

    if isa(W,'sdpvar')
        PSD_constraints = W >= 0; 
    else
        PSD_constraints = all(min(eig(W)) >= -tol);
    end
    
end