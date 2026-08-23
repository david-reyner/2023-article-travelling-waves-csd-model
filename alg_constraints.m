function [Vn, Va] = alg_constraints(Vn0, Va0, K_e, rhoN, rhoA, hp, conc)

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
    P_K = 1e-6; % astrocyte K+ permeability coefficient (cm/s)
    P_Na = 0.015e-6; % astrocyte Na+ permeability coefficient (cm/s)

    % Total initial sodium and potassium concentrations
    Natot = omega_n*Na_i0 + omega_a*Na_iA + omega_e*Na_e;
    Ktot = omega_n*K_i0 + omega_a*K_iA + omega_e*K_e0;
    
    % Intracellular sodium/potassium concentrations
    K_i = 1/omega_n*(Ktot - omega_e*K_e - omega_a*K_iA); % all(K_i == K_i0 + alpha0*(K_e0 - K_e))
    Na_i = 1/omega_n*(Natot - omega_e*Na_e - omega_a*Na_iA);
    
    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    E_Na = E(Na_i,Na_e); E_K = E(K_i,K_e);
    
    % Goldman-Hodgkin-Katz formulation
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    G1 = @(Phi,Ii,Ie) ((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    
    % Asymptotic functions
    mInf = @(v) 1./(1 + exp(-(v + 34)/5));
    nInf = @(v) 1./(1 + exp(-(v + 55)/14));
    mpInf = @(v) 1./(1 + exp(-(v + 40)/6));
    
    % Derivatives of asymptotic functions
    dmInf = @(v) (1/5)*exp(-(v + 34)/5)./(1 + exp(-(v + 34)/5)).^2;
    dnInf = @(v) (1/14)*exp(-(v + 55)/14)./(1 + exp(-(v + 55)/14)).^2;
    dmpInf = @(v) (1/6)*exp(-(v + 40)/6)./(1 + exp(-(v + 40)/6)).^2;

    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_Na = 7.7; % half activation Na+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes % <--------------
    
    % Pump formula (neurons)
    IPn = rhoN*(K_e./(K_K + K_e)).^2.*(Na_i./(K_Na + Na_i)).^3;
    
    % Pump formula (astrocytes)
    IPa = rhoA*(K_e./(K_K + K_e)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;
    
    % Ionic currents to neurons
    INa = @(v) gNa*mInf(v).^3.*(1 - nInf(v)).*(v - E_Na); % Fast sodium current
    INaP = @(v) gNaP*mpInf(v).*hp.*(v - E_Na); % Persistent sodium current
    IK = @(v) gK*(nInf(v).^4).*(v - E_K); % Persistent potassium current
    IL = @(v) gL*(v - E_L); % Leak current
    
    % Ionic currents to astrocytes
    phi_A = @(v) v*F/(R*T);
    INa_A = @(v) G(P_Na,phi_A(v),Na_iA,Na_e).*(abs(phi_A(v)) > 1e-10) + ...
                 P_Na*F*(Na_iA - Na_e).*(abs(phi_A(v)) <= 1e-10); % Fast sodium current
    IK_A = @(v) G(P_K,phi_A(v),K_iA,K_e).*(abs(phi_A(v)) > 1e-10) + ...
                P_K*F*(K_iA - K_e).*(abs(phi_A(v)) <= 1e-10); % Potassium current
    
    % Algebraic constraints
    cons1 = @(v) (INa(v) + INaP(v) + IK(v) + IL(v) + IPn); % Phi_Vn
    cons2 = @(v) (INa_A(v) + IK_A(v) + IPa); % Phi_Va
    
    % Derivatives of the algebraic constraints (wrt Vn and Va, resp.)
    dcons1 = @(v) gNa*(3*mInf(v).^2).*dmInf(v).*(1-nInf(v)).*(v - E_Na) + ...
        gNa*mInf(v).^3.*(1 - dnInf(v).*v - nInf(v) + dnInf(v)*E_Na) + ...
        gNaP*hp*dmpInf(v).*(v - E_Na) + gNaP*hp*mpInf(v) + ...
        gK*(4*nInf(v).^3).*dnInf(v).*(v - E_K) + gK*nInf(v).^4 + gL;
    dcons2 = @(v) P_Na*F*((1/phi)*G1(phi_A(v),Na_iA,Na_e) + ...
        phi_A(v).*((1/phi)*exp(-phi_A(v)).*(Na_e - Na_iA)./(exp(-phi_A(v)) - 1).^2)) + ...
        P_K*F*((1/phi)*G1(phi_A(v),K_iA,K_e) + ...
        phi_A(v).*((1/phi)*exp(-phi_A(v)).*(K_e - K_iA)./(exp(-phi_A(v)) - 1).^2));
    
    % Apply a Newton method to compute Vn and Va as functions of K_e
    [Vn, iter1, err1] = newton(Vn0, cons1, dcons1);
    [Va, iter2, err2] = newton(Va0, cons2, dcons2);
end

%_________________________________________________________________________%
function [x0, iter, err] = newton(x0, f, fp)
    % Newton method to compute, given an initial vector of seeds x0, the
    % zeros for the vectorial function f
    err = 1; % initial error
    iter = 1; % number of iterations
    while err >= 1e-12 && iter < 40
        x1 = x0 - f(x0)./fp(x0); % Newton correction
        err = norm(x1 - x0); % update error
        x0 = x1; % update seed
        iter = iter + 1;
    end
end