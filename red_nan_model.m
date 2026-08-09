function dx = red_nan_model(t, x, N, rhoN, rhoA, hp, conc)

    global pumpK

    % Fixed sodium and potassium concentrations
    Na_e = conc(1); Na_i0 = conc(2); Na_iA = conc(3);
    K_e0 = conc(4); K_i0 = conc(5); K_iA = conc(6);
    
    % Neurons' membrane potential
    Vn = x(1:N);
    
    % Extracellular concentration of potassium
    K_e = x(N+1:2*N);
    
    % Astrocytes' membrane potential
    Va = x(2*N+1:3*N);
    
    % Neurons and astrocytes parameters
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    F = 96485; % Faraday constant (coulomb/M)
    T = 310; % absolute temperature (K)
    phi = R*T/F;
    Cm = 1; % neuron membrane capacitance per unit area (microF/cm2)
    Cma = 1; % astrocyte membrane capacitance per unit area (microF/cm2)
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
    omega_e = alpha0*(omega_n + omega_a); % extracellular volume
    DK = 0.001; % K+ diffusion coefficient (scaled by square of neuronal distance)
    P_K = 1e-6; % astrocyte K+ permeability coefficient (cm/s)
    P_Na = 0.015e-6; % astrocyte Na+ permeability coefficient (cm/s)

    % Total initial sodium and potassium concentrations
    Natot = omega_n*Na_i0 + omega_a*Na_iA + omega_e*Na_e;
    Ktot = omega_n*K_i0 + omega_a*K_iA + omega_e*K_e0;
    
    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    
    % Goldman-Hodgkin-Katz formulation
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    
    % Asymptotic functions
    mInf = 1./(1 + exp(-(Vn + 34)/5));
    nInf = 1./(1 + exp(-(Vn + 55)/14)); n = nInf; % instantaneous convergence to nInf
    mpInf = 1./(1 + exp(-(Vn + 40)/6));
    
    % Intracellular sodium/potassium concentrations
    K_i = 1/omega_n*(Ktot - omega_e*K_e - omega_a*K_iA); % all(K_i == K_i0 + alpha0*(K_e0 - K_e))
    Na_i = 1/omega_n*(Natot - omega_e*Na_e - omega_a*Na_iA);
    
    K_K = 2; % half activation K+ concentration (mM) for neurons
    K_Na = 7.7; % half activation Na+ concentration (mM) for neurons
    K_NaA = 7.7; % half activation Na+ concentration (mM) for astrocytes
    
    % Pump formula (neurons)
    IPn = rhoN*(K_e./(K_K + K_e)).^2.*(Na_i./(K_Na + Na_i)).^3;
    
    % Pump formula (astrocytes)
    IPa = rhoA*(K_e./(K_K + K_e)).^2.*(Na_iA./(K_NaA + Na_iA)).^3;

    % Ionic currents to neurons
    INa = gNa*mInf.^3.*(1 - n).*(Vn - E(Na_i,Na_e)); % Fast sodium current
    INaP = gNaP*mpInf.*hp.*(Vn - E(Na_i,Na_e)); % Persistent sodium current
    IK = gK*(n.^4).*(Vn - E(K_i,K_e)); % Persistent potassium current
    IL = gL*(Vn - E_L); % Leak current
    
    % Ionic currents to astrocytes
    IK_A = zeros(N,1); INa_A = zeros(N,1);
    for i = 1:N
        phi_A = Va(i)*F/(R*T);
        if abs(phi_A) > 1e-10
            INa_A(i) = G(P_Na,phi_A,Na_iA,Na_e); % Fast sodium current
            IK_A(i) = G(P_K,phi_A,K_iA,K_e(i)); % Potassium current
        else
            INa_A(i) = P_Na*F*(Na_iA - Na_e); % Fast sodium current
            IK_A(i) = P_K*F*(K_iA - K_e(i)); % Potassium current
        end
    end
    
    % Dirichlet boundary conditions on extracellular potassium around 
    % first/last neurons
    K_eN = K_e0;
    
    % Injection of K+ to the extracellular space of the middle cells
    Iinj = zeros(N,1);
    mid = [N/2-1 N/2 N/2+1 N/2+2]; % 24,25,26,27
    if pumpK
        Iinj(mid) = 0.005;
    end
    
    % Reduced model for a neuron-astrocyte network responsible for 
    % initiation of cortical spreading depolarization
    dx = zeros(3*N,1);
    dx(1:N) = -(1/Cm)*(INa + INaP + IK + IL + IPn);

    dx(N+1) = DK*(K_e(2) - 2*K_e(1) + K_e0) ...
              + 10*Sn/(F*omega_e)*(IK(1) - 2*IPn(1)) ...
              + 10*Sa/(F*omega_e)*(IK_A(1) - 2*IPa(1));
    dx(N+2:2*N-1) = DK*(K_e(3:N) - 2*K_e(2:N-1) + K_e(1:N-2)) ...
                    + 10*Sn/(F*omega_e)*(IK(2:N-1) - 2*IPn(2:N-1)) + Iinj(2:N-1) ...
                    + 10*Sa/(F*omega_e)*(IK_A(2:N-1) - 2*IPa(2:N-1));
    dx(2*N) = DK*(K_e(N-1) - 2*K_e(N) + K_eN) ...
              + 10*Sn/(F*omega_e)*(IK(N) - 2*IPn(N)) ...
              + 10*Sa/(F*omega_e)*(IK_A(N) - 2*IPa(N));
    
    dx(2*N+1:3*N) = -(1/Cma)*(INa_A + IK_A + IPa);
    
    if pumpK
        if all(Vn(mid) > -30)
            t
            pumpK = false;
        end
    end
end