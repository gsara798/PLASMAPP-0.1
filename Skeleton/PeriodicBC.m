%% Function for periodic boundary conditions
function variable = PeriodicBC(variable,LB,UB)
    % LB = Lower boundary, UB = Upper Boundary
    out=(variable<LB); variable(out) = variable(out) + UB;
    out=(variable>UB); variable(out) = variable(out) - UB;
end
