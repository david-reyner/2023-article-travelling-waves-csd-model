function [varargout] = diff_PhiVa(x0, rhoA, conc, varargin)
    
    % Transpose if x0 is a column vector
    if size(x0,2) == 1
        x0 = x0';
    end
    
    % Initial point
    Va = x0(:,1); K_e = x0(:,2);
    
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
    
    % Goldman-Hodgkin-Katz formulation (old implementation)
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    G1 = @(Phi,Ii,Ie) ((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes
                
    % Pump formula (astrocytes)
    IPa = rhoA*(K_e./(K_K + K_e)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;

    % Ionic currents to astrocytes
    phi_A = Va*F/(R*T);
%     if abs(phi_A) > 1e-10
%         INa_A = G(P_Na,phi_A,Na_iA,Na_e); % Fast sodium current
%         IK_A = G(P_K,phi_A,K_iA,K_e); % Potassium current
%     else
%         INa_A = P_Na*F*(Na_iA - Na_e); % Fast sodium current
%         IK_A = P_K*F*(K_iA - K_e); % Potassium current
%     end
    INa_A = GHK_formula(P_Na, F, phi_A, Na_iA, Na_e*ones(size(phi_A))); % Fast sodium current
    IK_A = GHK_formula(P_K, F, phi_A, K_iA, K_e); % Potassium current
    
    % Indices where definition for the GHK formula applies
    id = abs(phi_A) > 1e-10;

    % Phi_Va
    Phi_Va = -(INa_A + IK_A + IPa);
        
    % First derivatives of Phi_Va with respect to Va and [K^+]_e
    dVPhi_Va = -(P_Na*F*((1/phi)*G1(phi_A,Na_iA,Na_e) + ...
        phi_A.*((1/phi)*exp(-phi_A).*(Na_e - Na_iA)./(exp(-phi_A) - 1).^2)) + ...
        P_K*F*((1/phi)*G1(phi_A,K_iA,K_e) + ...
        phi_A.*((1/phi)*exp(-phi_A).*(K_e - K_iA)./(exp(-phi_A) - 1).^2)));
    dVPhi_Va(~id) = -0.5*(F^2/(R*T))*(P_Na*(Na_iA + Na_e) + P_K*(K_iA + K_e(~id)));
    
    dKPhi_Va = -(P_K*F*phi_A.*(exp(-phi_A)./(exp(-phi_A) - 1)) + ...
        (2*rhoA*K_K*K_e./(K_e + K_K).^3).*(Na_iA./(K_NaA + Na_iA)).^3);
    dKPhi_Va(~id) = -(-P_K*F + ...
        (2*rhoA*K_K*K_e(~id)./(K_e(~id) + K_K).^3).*(Na_iA./(K_NaA + Na_iA)).^3);
    
    % Output
    if nargin > 3
        % When required, return only Phi_Va, dVPhi_Va or dKPhi_Va
        var = varargin{1};
        if var == 1
            varargout{1} = Phi_Va;
        elseif var == 2
            varargout{1} = dVPhi_Va;
        else
            varargout{1} = dKPhi_Va;
        end
    else
        % Return the function Phi_Va and its derivatives
        varargout{1} = Phi_Va;
        varargout{2} = dVPhi_Va;
        varargout{3} = dKPhi_Va;
    end
end

%_________________________________________________________________________%
function I = GHK_formula(P, F, Phi, Ii, Ie)
    % Goldman-Hodgkin-Katz formulation in vectorized form
    id = abs(Phi) > 1e-10;
    I = zeros(size(Phi));
    I(id) = P*F*Phi(id).*((Ie(id).*exp(-Phi(id)) - Ii)./(exp(-Phi(id)) - 1)); % abs(Phi) > 1e-10
    I(~id) = P*F*(Ii - Ie(~id)); % abs(Phi) <= 1e-10
end