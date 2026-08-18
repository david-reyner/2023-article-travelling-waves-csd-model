function Vn = Phi_Vn(Vn0, K_e, rhoN, hp, conc)

    % Function Phi_Vn
    f = @(v) diff_PhiVn([v K_e], rhoN, hp, conc, 1);
    
    % Derivative of function Phi_Vn with respect to Vn
    df = @(v) diff_PhiVn([v K_e], rhoN, hp, conc, 2);
    
    % Apply a Newton method to compute Vn(K_e) from Phi_Vn = 0
    [Vn, iter, err] = newton(Vn0, f, df);
end

%_________________________________________________________________________%
function [x0, iter, err] = newton(x0, f, fp)
    % Vectorized Newton method to compute, given an initial vector of seeds
    % x0, the zeros of a scalar function f
    err = ones(length(x0),1); % initial error
    iter = zeros(length(x0),1); % number of iterations
    id = err >= 1e-12 & iter < 40; % indices of seeds that have not yet converged
    while any(id) % err(id) >= 1e-12 && iter(id) < 40
        delta_x = -f(x0)./fp(x0);
        x1 = x0(id) + delta_x(id); % Newton correction
        err(id) = abs(x1 - x0(id)); % update error
        x0(id) = x1; % update seed
        iter(id) = iter(id) + 1; % update num. iterations
        id = err >= 1e-12 & iter < 40; % update indices
    end
end