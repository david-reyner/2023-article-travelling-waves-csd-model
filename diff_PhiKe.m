function [varargout] = diff_PhiKe(x0, rhoN, rhoA, conc, varargin)
    
    % Transpose if x0 is a column vector
    if size(x0,2) == 1
        x0 = x0';
    end

    % Initial point
    Vn = x0(:,1); Va = x0(:,2); K_e = x0(:,3);

    % Fixed sodium and potassium concentrations
    Na_e = conc(1); Na_i0 = conc(2); Na_iA = conc(3);
    K_e0 = conc(4); K_i0 = conc(5); K_iA = conc(6);
    
    % Neurons and astrocytes parameters
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    F = 96485; % Faraday constant (coulomb/M)
    T = 310; % absolute temperature (K)
    phi = R*T/F;
    gK = 15; % maximum delayed rectifier potassium current conductance (mS/cm2)
    Sn = 922; % neuron total membrane surface area (micro m2)
    Sa = 1600; % neuron total membrane surface area (micro m2)
    omega_n = 2160; % neuron total volume
    omega_a = 2000; % astrocyte total volume
    alpha0 = 0.2; % volume fraction between neuron/astrocyte and extracel. space
    omega_e = alpha0*(omega_n + omega_a); % extracellular volume
    beta = omega_a/omega_n; % ratio between neuron and astrocyte volumes
    P_K = 1e-6; % astrocyte K+ permeability coefficient (cm/s)
    
    % Total initial sodium and potassium concentrations
    Natot = omega_n*Na_i0 + omega_a*Na_iA + omega_e*Na_e;
    Ktot = omega_n*K_i0 + omega_a*K_iA + omega_e*K_e0;
    
    % Intracellular sodium/potassium concentrations
    K_i = 1/omega_n*(Ktot - omega_e*K_e - omega_a*K_iA); % all(K_i == K_i0 + alpha0*(K_e0 - K_e))
    Na_i = 1/omega_n*(Natot - omega_e*Na_e - omega_a*Na_iA);
    
    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    E_K = E(K_i,K_e);
    
    % Goldman-Hodgkin-Katz formulation
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    G1 = @(Phi,Ii,Ie) ((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));

    % n-asymptotic function
    nInf = 1./(1 + exp(-(Vn + 55)/14));
    
    % Derivatives of the n asymptotic function
    dnInf = (1/14)*exp(-(Vn + 55)/14)./(1 + exp(-(Vn + 55)/14)).^2;

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_Na = 7.7; % half activation Na+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes
    
    % Pump formula (neurons)
    IPn = rhoN*(K_e./(K_K + K_e)).^2.*(Na_i./(K_Na + Na_i)).^3;
    
    % Pump formula (astrocytes)
    IPa = rhoA*(K_e./(K_K + K_e)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;
    
    % Persistent potassium current to neurons
    IK = gK*(nInf.^4).*(Vn - E_K);
    
    % Potassium current to astrocytes
    phi_A = Va*F/(R*T);
    IK_A = GHK_formula(P_K, F, phi_A, K_iA, K_e); % Potassium current

    % Indices where different definition of the GHK formula applies
    id = abs(phi_A) > 1e-10;
    
    % Phi_Ke
    Phi_Ke = 10*Sn/(F*omega_e)*(IK - 2*IPn) + ...
             10*Sa/(F*omega_e)*(IK_A - 2*IPa);

    % First derivatives of the Phi_Ke with respect to Vn, Va and [K^+]_e
    dVnPhi_Ke = 10*Sn/(F*omega_e)*(gK*(4*nInf.^3).*dnInf.*(Vn - E_K) ...
                + gK*nInf.^4);
    
    dVaPhi_Ke = 10*Sa/(F*omega_e)*(P_K*F*((1/phi)*G1(phi_A,K_iA,K_e) + ...
                    phi_A.*((1/phi)*exp(-phi_A).*(K_e - K_iA)./(exp(-phi_A) - 1).^2)));
    dVaPhi_Ke(~id) = 0;     
    
    dKePhi_Ke = 10*Sn/(F*omega_e)*(-gK*nInf.^4.*phi.*(K_i + alpha0*(1 + beta)*K_e)./(K_e.*K_i) - ...
                2*(2*rhoN*K_K*K_e./(K_e + K_K).^3).*(Na_i./(K_Na + Na_i)).^3) + ...
                10*Sa/(F*omega_e)*(P_K*F*phi_A.*(exp(-phi_A)./(exp(-phi_A) - 1)) - ...
                2*(2*rhoA*K_K*K_e./(K_e + K_K).^3).*(Na_iA./(K_NaA + Na_iA)).^3);
    dKePhi_Ke(~id) = 10*Sn/(F*omega_e)*(-gK*nInf(~id).^4.*phi.*(K_i(~id) + alpha0*(1 + beta)*K_e(~id))./(K_e(~id).*K_i(~id)) - ...
                     2*(2*rhoN*K_K*K_e(~id)./(K_e(~id) + K_K).^3).*(Na_i./(K_Na + Na_i)).^3) + ...
                     10*Sa/(F*omega_e)*(-P_K*F + ...
                     2*(2*rhoA*K_K*K_e(~id)./(K_e(~id) + K_K).^3).*(Na_iA./(K_NaA + Na_iA)).^3);
    
    % Output
    if nargin > 4
        % When required, return only Phi_Ke, dVnPhi_Ke, dVaPhi_Ke or dKePhi_Ke
        var = varargin{1};
        if var == 1
            varargout{1} = Phi_Ke;
        elseif var == 2
            varargout{1} = dVnPhi_Ke;
        elseif var == 3
            varargout{1} = dVaPhi_Ke;
        else
            varargout{1} = dKePhi_Ke;
        end
    else
        % Return the function Phi_Va and its derivatives
        varargout{1} = Phi_Ke;
        varargout{2} = dVnPhi_Ke;
        varargout{3} = dVaPhi_Ke;
        varargout{4} = dKePhi_Ke;
    end
end

%_________________________________________________________________________%
function I = GHK_formula(P, F, Phi, Ii, Ie)
    % Goldman-Hodgkin-Katz formulation (vectorized)
    id = abs(Phi) > 1e-10;
    I = zeros(size(Phi));
    I(id) = P*F*Phi(id).*((Ie(id).*exp(-Phi(id)) - Ii)./(exp(-Phi(id)) - 1)); % abs(Phi) > 1e-10
    I(~id) = P*F*(Ii - Ie(~id)); % abs(Phi) <= 1e-10
end