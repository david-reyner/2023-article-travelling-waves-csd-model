function Ke0c = cont_PhiVn2(x0, sgn, rhoN, hp, conc)
    % Implementation of the Keller continuation method
    
    % x0: initial solution on the curve to be continuated
    % xd: initial direction at solution x0
    
    % Compute an initial direction at the initial solution x0
    [~, dVPhi_Vn, dKPhi_Vn] = diff_PhiVn(flip(x0), rhoN, hp, conc);
    xd = [-dVPhi_Vn/dKPhi_Vn; 1]/norm([-dVPhi_Vn/dKPhi_Vn; 1]);
    
    % Zero curve of Phi_Vn: Ke in terms of Vn
    Ke0c = [x0' 0 0];
    
    % Step on the tangent lines (sign determines direction to follow)
    delta_s = sgn*0.02;
    
    while (sgn < 0 && x0(2) > -80) || (sgn > 0 && x0(1) < 311)
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Newton's method %%%%%%%%%%%%%%%%%%%%%%%
        % Prediction
        xaux = x0 + delta_s*xd;
        
        % Correction
        error = 1; % initial error
        iter = 1; % number of iterations
        while error >= 1e-12 && iter < 20
            [Phi_Vn, dVPhi_Vn, dKPhi_Vn] = diff_PhiVn(flip(xaux), rhoN, hp, conc);
            A = [dKPhi_Vn dVPhi_Vn;
                   xd(1)    xd(2)];
            b = [Phi_Vn; (xaux - x0)'*xd - delta_s];
            delta_x = -A\b;
            x1 = xaux + delta_x; % Newton correction
            error = norm(b); % error in the target function
            errStep = norm(x1 - xaux); % update error on the step
            xaux = x1; % update approximation
            iter = iter + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if error < 1e-12
            % Restore initial step size
            delta_s = sgn*0.02;
            
            % Update initial direction xd
            [~, dVPhi_Vn, dKPhi_Vn] = diff_PhiVn(flip(xaux), rhoN, hp, conc);
            A = [dKPhi_Vn dVPhi_Vn;
                   xd(1)    xd(2)];
            xd = A\[0; 1]; % tangent vector to the curve preserving orientation

            % Update initial solution x0
            x0 = xaux;

            Ke0c = [Ke0c; x0' error iter];
        elseif abs(delta_s) > 0.001
            % Reduced step size along tangent line
            fprintf("\nStep size reduced from %1.5f to %1.5f\n", delta_s, delta_s/2);
            delta_s = delta_s/2;
        else % Stop continuation method
            x0 = [500, -85];
            disp("Step along tangent line too small.")
        end
    end
end