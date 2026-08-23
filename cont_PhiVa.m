function Va0c = cont_PhiVa(x0, rhoA, conc)
    % Implementation of a simple continuation method
        
    % Fixed sodium and potassium concentrations
    Na_e = conc(1); Na_i0 = conc(2); Na_iA = conc(3);
    K_e0 = conc(4); K_i0 = conc(5); K_iA = conc(6);
    
    % Parameters
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    F = 96485; % Faraday constant (coulomb/M)
    T = 310; % absolute temperature (K)
    phi = R*T/F;
    P_K = 1e-6; % astrocyte K+ permeability coefficient (cm/s)
    P_Na = 0.015e-6; % astrocyte Na+ permeability coefficient (cm/s)
    
    % Goldman-Hodgkin-Katz formulation
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    G1 = @(Phi,Ii,Ie) ((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes

    % Parameter of the curve
    Ke = linspace(0, 600, 6500);

    % Zero curve of Phi_Va
    Va0c = zeros(length(Ke), 4);

    for i = 1:length(Ke)
        
        K_e = Ke(i);
        
        % Pump formula (astrocytes)
        IPa = rhoA*(K_e./(K_K + K_e)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;
        
        % Ionic currents to astrocytes
        INa_A = @(v) GHK_formula(v, P_Na, F, Na_iA, Na_e); % Fast sodium current
        IK_A = @(v) GHK_formula(v, P_K, F, K_iA, K_e); % Potassium current
%         INa_A = @(v) G(P_Na,phi_A(v),Na_iA,Na_e).*(abs(phi_A(v)) > 1e-10) + ...
%                  P_Na*F*(Na_iA - Na_e).*(abs(phi_A(v)) <= 1e-10); 
%         IK_A = @(v) G(P_K,phi_A(v),K_iA,K_e).*(abs(phi_A(v)) > 1e-10) + ...
%                 P_K*F*(K_iA - K_e).*(abs(phi_A(v)) <= 1e-10); 
        
        % Phi_Va
        Phi_Va = @(v) (INa_A(v) + IK_A(v) + IPa);
        
        % Derivative of Phi_Va with respect to Va
        phi_A = @(v) v*F/(R*T);
        dPhi_Va = @(v) P_Na*F*((1/phi)*G1(phi_A(v),Na_iA,Na_e) + ...
            phi_A(v).*((1/phi)*exp(-phi_A(v)).*(Na_e - Na_iA)./(exp(-phi_A(v)) - 1).^2)) + ...
            P_K*F*((1/phi)*G1(phi_A(v),K_iA,K_e) + ...
            phi_A(v).*((1/phi)*exp(-phi_A(v)).*(K_e - K_iA)./(exp(-phi_A(v)) - 1).^2));
        
        % Newton's method
        error = 1; % initial error
        iter = 1; % number of iterations
        while error >= 1e-12 && iter < 60
            x1 = x0 - Phi_Va(x0)./dPhi_Va(x0); % Newton correction
            error = norm(x1 - x0); % update error
            x0 = x1; % update approximation
            iter = iter + 1;
        end
        
        Va0c(i,:) = [x0 K_e error iter]; 
    end
end

%_________________________________________________________________________%
function I = GHK_formula(v, P, F, Ii, Ie)
    % Constants
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    T = 310; % absolute temperature (K)

    % Goldman-Hodgkin-Katz formulation (vectorized)
    Phi = v*F/(R*T);
    id = abs(Phi) > 1e-10;
    I = zeros(size(Phi));
    I(id) = P*F*Phi(id).*((Ie(id).*exp(-Phi(id)) - Ii)./(exp(-Phi(id)) - 1)); % abs(Phi) > 1e-10
    I(~id) = P*F*(Ii - Ie(~id)); % abs(Phi) <= 1e-10
end