function fun = rhs_red_nan(x, y, z, rhoN, rhoA, hp, conc, varargin)
    % Return the right hand side of the reduced neuron-astrocyte network
    % model, that is, functions Phi_Vn, Phi_Va or Phi_Ke
    
    % Treat x and y as approximations of the curves x = x(z) and y = y(z)
    % to be refined via a Newton method
    if nargin < 9
        x = Phi_Vn(x, z, rhoN, hp, conc);
        y = Phi_Va(y, z, 0, rhoA, conc);
    end

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
    Sn = 922; % neuron total membrane surface area (micro m2)
    Sa = 1600; % neuron total membrane surface area (micro m2)
    omega_n = 2160; % neuron total volume
    omega_a = 2000; % astrocyte total volume
    alpha0 = 0.2; % volume fraction between neuron/astrocyte and extracel. space
    omega_e = alpha0*(omega_n + omega_a); % extracelular volume
    P_K = 1e-6; % astrocyte K+ permeability coefficient (cm/s)
    P_Na = 0.015e-6; % astrocyte Na+ permeability coefficient (cm/s)

    % Total initial sodium and potassium concentrations
    Natot = omega_n*Na_i0 + omega_a*Na_iA + omega_e*Na_e;
    Ktot = omega_n*K_i0 + omega_a*K_iA + omega_e*K_e0;
    
    % Intracellular sodium/potassium concentrations
    K_i = 1/omega_n*(Ktot - omega_e*z - omega_a*K_iA); % all(K_i == K_i0 + alpha0*(K_e0 - K_e))
    Na_i = 1/omega_n*(Natot - omega_e*Na_e - omega_a*Na_iA);
    
    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    E_Na = E(Na_i,Na_e); E_K = E(K_i,z);
        
    % Asymptotic functions
    mInf = 1./(1 + exp(-(x + 34)/5));
    nInf = 1./(1 + exp(-(x + 55)/14));
    mpInf = 1./(1 + exp(-(x + 40)/6));

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_Na = 7.7; % half activation Na+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes
    
    % Pump formula (neurons)
    IPn = rhoN*(z./(K_K + z)).^2.*(Na_i./(K_Na + Na_i)).^3;
    
    % Pump formula (astrocytes)
    IPa = rhoA*(z./(K_K + z)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;
    
    % Ionic currents to neurons
    INa = gNa*mInf.^3.*(1 - nInf).*(x - E_Na); % Fast sodium current
    INaP = gNaP*mpInf.*hp.*(x - E_Na); % Persistent sodium current
    IK = gK*(nInf.^4).*(x - E_K); % Persistent potassium current
    IL = gL*(x - E_L); % Leak current
    
    % Ionic currents to astrocytes
    phi_A = y*F/(R*T);
    INa_A = GHK_formula(P_Na, F, phi_A, Na_iA, Na_e*ones(size(phi_A))); % Fast sodium current
    IK_A = GHK_formula(P_K, F, phi_A, K_iA, z); % Potassium current
%     INa_A = @(v) G(P_Na,phi_A(v),Na_iA,Na_e).*(abs(phi_A(v)) > 1e-10) + ...
%                  P_Na*F*(Na_iA - Na_e).*(abs(phi_A(v)) <= 1e-10); % Fast sodium current
%     IK_A = @(v) G(P_K,phi_A(v),K_iA,z).*(abs(phi_A(v)) > 1e-10) + ...
%                 P_K*F*(K_iA - z).*(abs(phi_A(v)) <= 1e-10); % Potassium current
    
    % Return functions Phi_Vn, Phi_Va and Phi_Ke (when required)
    var = varargin{1};
    if var == 1
        fun = -(INa + INaP + IK + IL + IPn); % Phi_Vn
    elseif var == 2
        fun = -(INa_A + IK_A + IPa); % Phi_Va
    else
        fun = 10*Sn/(F*omega_e)*(IK - 2*IPn) + ...
              10*Sa/(F*omega_e)*(IK_A - 2*IPa); % Phi_Ke
    end
    
    if ~isreal(fun)
        fun = real(fun);
    end
end

%_________________________________________________________________________%
function I = GHK_formula(P, F, Phi, Ii, Ie)
    % Goldman-Hodgkin-Katz formulation in vectorized form
%     G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    id = abs(Phi) > 1e-10;
    I = zeros(size(Phi));
    I(id) = P*F*Phi(id).*((Ie(id).*exp(-Phi(id)) - Ii)./(exp(-Phi(id)) - 1)); % abs(Phi) > 1e-10
    I(~id) = P*F*(Ii - Ie(~id)); % abs(Phi) <= 1e-10
end