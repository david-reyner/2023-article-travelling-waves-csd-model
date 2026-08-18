function Va = Phi_Va(Va0, K_e, m, rhoA, conc)

    % Function Phi_Va = m (m constant)
    g = @(v) diff_PhiVa([v K_e], rhoA, conc, 1) - m;
    
    % Derivative of function Phi_Va - m with respect to Va
    dg = @(v) diff_PhiVa([v K_e], rhoA, conc, 2);
    
    % Apply a Newton method to compute Va(K_e) from Phi_Va = m (m ~= 0)
    [Va, iter, err] = newton(Va0, g, dg);
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