function [x1, varargout] = secant_method(F, x0, x1, tol)
    % Secant method to find a zero of F(x) with initial values x0 and x1 
    % with a given tolerance tol
    iter = 1; % Current iteration
    maxiter = 20; % Max. number of iterations
    while iter <= maxiter && abs(x1 - x0) >= tol 
%         abs(x1 - x0)
        f0 = F(x0); f1 = F(x1); % evaluating function F at points x0 and x1
        x2 = x1 - f1*(x1 - x0)/(f1 - f0); % approximation by secant line
        
        % Update initial values
        x0 = x1; x1 = x2;
        
        % Increment iteration counter
        iter = iter + 1;
    end
    
    % If more than one output is required, return error and number of
    % iterations
    if nargout > 1
        varargout{1} = abs(x1 - x0);
        varargout{2} = iter;
    end
end