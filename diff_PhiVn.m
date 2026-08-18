function [varargout] = diff_PhiVn(x0, rhoN, hp, conc, varargin)
    
    % Transpose if x0 is a column vector
    if size(x0,2) == 1
        x0 = x0';
    end

    % Initial point
    Vn = x0(:,1); K_e = x0(:,2);

    % Fixed sodium and potassium concentrations
    Na_e = conc(1); Na_i0 = conc(2); Na_iA = conc(3);
    K_e0 = conc(4); K_i0 = conc(5); K_iA = conc(6);
    
    % Neurons and astrocytes parameters
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    F = 96485; % Faraday constant (coulomb/M)
    T = 310; % absolute temperature (K)
    phi = R*T/F;
    gNa = 50.0; % maximum fast sodium current conductance (mS/cm2)
    gNaP = 0.8; % maximum persistent sodium current conductance (mS/cm2)
    gK = 15; % maximum delayed rectifier potassium current conductance (mS/cm2)
    gL = 0.5; % leak current conductance (mS/cm2)
    E_L = -70; % leak current reversal potential (mV)
    omega_n = 2160; % neuron total volume
    omega_a = 2000; % astrocyte total volume
    alpha0 = 0.2; % volume fraction between neuron/astrocyte and extracel. space
    omega_e = alpha0*(omega_n + omega_a); % extracellular volume
    beta = omega_a/omega_n; % ratio between neuron and astrocyte volumes

    % Total initial sodium and potassium concentrations
    Natot = omega_n*Na_i0 + omega_a*Na_iA + omega_e*Na_e;
    Ktot = omega_n*K_i0 + omega_a*K_iA + omega_e*K_e0;
    
    % Intracellular sodium/potassium concentrations
    K_i = 1/omega_n*(Ktot - omega_e*K_e - omega_a*K_iA); % all(K_i == K_i0 + alpha0*(K_e0 - K_e))
    Na_i = 1/omega_n*(Natot - omega_e*Na_e - omega_a*Na_iA);
    
    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    E_Na = E(Na_i,Na_e); E_K = E(K_i,K_e);

    % Asymptotic functions
    mInf = 1./(1 + exp(-(Vn + 34)/5));
    nInf = 1./(1 + exp(-(Vn + 55)/14));
    mpInf = 1./(1 + exp(-(Vn + 40)/6));
    
    % Derivatives of asymptotic functions
    dmInf = (1/5)*exp(-(Vn + 34)/5)./(1 + exp(-(Vn + 34)/5)).^2;
    dnInf = (1/14)*exp(-(Vn + 55)/14)./(1 + exp(-(Vn + 55)/14)).^2;
    dmpInf = (1/6)*exp(-(Vn + 40)/6)./(1 + exp(-(Vn + 40)/6)).^2;

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_Na = 7.7; % half activation Na+ concentration (mM) for neurons
    
    % Pump formula (neurons)
    IPn = rhoN*(K_e./(K_K + K_e)).^2.*(Na_i./(K_Na + Na_i)).^3;
    
    % Ionic currents to neurons
    INa = gNa*mInf.^3.*(1 - nInf).*(Vn - E_Na); % Fast sodium current
    INaP = gNaP*mpInf.*hp.*(Vn - E_Na); % Persistent sodium current
    IK = gK*(nInf.^4).*(Vn - E_K); % Persistent potassium current
    IL = gL*(Vn - E_L); % Leak current

    % Phi_Vn
    Phi_Vn = -(INa + INaP + IK + IL + IPn);
    
    % First derivatives of Phi_Vn with respect to Vn and [K^+]_e
    dVPhi_Vn = -(gNa*(3*mInf.^2).*dmInf.*(1-nInf).*(Vn - E_Na) + ...
        gNa*mInf.^3.*(1 - dnInf.*Vn - nInf + dnInf*E_Na) + ...
        gNaP*hp*dmpInf.*(Vn - E_Na) + gNaP*hp*mpInf + ...
        gK*(4*nInf.^3).*dnInf.*(Vn - E_K) + gK*nInf.^4 + gL);
    dKPhi_Vn = -(-gK*nInf.^4.*phi.*(K_i + alpha0*(1 + beta)*K_e)./(K_e.*K_i) + ...
        (2*rhoN*K_K*K_e./(K_e + K_K).^3).*(Na_i./(K_Na + Na_i)).^3);
    
    % Output
    if nargin > 4
        % When required, return only Phi_Vn, dVPhi_Vn or dKPhi_Vn
        var = varargin{1};
        if var == 1
            varargout{1} = Phi_Vn;
        elseif var == 2
            varargout{1} = dVPhi_Vn;
        else
            varargout{1} = dKPhi_Vn;
        end
    else
        % Return the function Phi_Vn and its derivatives
        varargout{1} = Phi_Vn;
        varargout{2} = dVPhi_Vn;
        varargout{3} = dKPhi_Vn;
    end
end