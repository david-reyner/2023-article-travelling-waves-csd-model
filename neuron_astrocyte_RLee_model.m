function dx = neuron_astrocyte_RLee_model(t, x, delta_x, N, eps, rhoN, rhoA, ...
                                          sig_gap, Ncn)
    
    global pumpK

    % Neurons' membrane potential
    Vn = x(1:N);

    % Gating variables n and hp
    n = x(N+1:2*N); hp = x(2*N+1:3*N);
    
    % Sodium and potassium concentrations in neurons' cytosol
    Na_e = x(3*N+1:4*N); Na_i = x(4*N+1:5*N);
    
    % Extracellular concentrations of sodium and potassium
    K_e = x(5*N+1:6*N); K_i = x(6*N+1:7*N);
    
    % Astrocytes' membrane potential
    Va = x(7*N+1:8*N);
    
    % Sodium and potassium concentrations in astrocytes' cytosol
    Na_iA = x(8*N+1:9*N); K_iA = x(9*N+1:10*N);
    
    % Neurons and astrocytes parameters
    R = 8310; % universal gas constant (J/M K) (factor 1000 to convert V to mV)
    F = 96485; % Faraday constant (coulomb/M)
    T = 310; % absolute temperature (K)
    phi = R*T/F;
    Cm = 1; % membrane capacitance per unit area (microF/cm2)
    Cma = 1; % membrane capacitance per unit area (microF/cm2)
    phi_n = 0.8; % maximum rate of activation of potassium channels (1/ms)
    phi_hp = 0.05; % maximum rate of activation of sodium channels (1/ms)
    gNa = 50.0; % maximum fast sodium current conductance (mS/cm2)
    gNaP = 0.8; % maximum persistent sodium current conductance (mS/cm2)
    gK = 15.0; % maximum delayed rectifier potassium current conductance (mS/cm2)
    gL = 0.5; % leak current conductance (mS/cm2)
    E_L = -70; % leak current reversal potential (mV)
    Sn = 922; % neuron total membrane surface area (micro m2)
    omega_n = 2160; % neuron total volume
    Sa = 1600; % astrocyte total membrane surface area (micro m2)
    omega_a = 2000; % astrocyte total volume
    alpha0 = 0.2; % volume fraction between neuron/astrocyte and extracel. space
    omega_e = alpha0*(omega_n + omega_a); % extracellular volume
    DNa = 0.000687; % Na+ diffusion coefficient (scaled by distance between neurons)
    DK = 0.001; % K+ diffusion coefficient (scaled by distance between neurons)
    P_K = 1e-6; % membrane K+ permeability coefficient (cm/s)
    P_Na = 0.015e-6; % membrane Na+ permeability coefficient (cm/s)

    % Reversal potentials given by the Nernst equation
    E = @(Ii,Ie) phi*log(Ie./Ii);
    
    % Goldman-Hodgkin-Katz formulation
    G = @(P,Phi,Ii,Ie) P*F*Phi.*((Ie.*exp(-Phi) - Ii)./(exp(-Phi) - 1));
    
    % Asymptotic functions
    mInf = 1./(1 + exp(-(Vn + 34)/5));
    nInf = 1./(1 + exp(-(Vn + 55)/14));
    mpInf = 1./(1 + exp(-(Vn + 40)/6));
    hpInf = 1./(1 + exp(-(Vn + 48)/(-6)));

    % Theta functions
    theta_n = 0.05 + 0.27./(1 + exp(-(Vn + 40)/(-12)));
    theta_hp = 10000./cosh((Vn + 48)/12);
    
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
            INa_A(i) = G(P_Na,phi_A,Na_iA(i),Na_e(i)); % Fast sodium current
            IK_A(i) = G(P_K,phi_A,K_iA(i),K_e(i)); % Potassium current
        else
            INa_A(i) = P_Na*F*(Na_iA(i) - Na_e(i)); % Fast sodium current
            IK_A(i) = P_K*F*(K_iA(i) - K_e(i)); % Potassium current
        end
    end
    
    % Currents due to gap junctions among astrocytes
    P_K_gap = sig_gap*P_K; % K+ gap junction permeability
    P_Na_gap = 0.8*P_K_gap; % Na+ gap junction permeability
    IK_gap = zeros(N,1); INa_gap = zeros(N,1);
    for j = 1:N
        for k = max(j-Ncn,1):min(j+Ncn,N) % loop on number of gap connections per side
            if k ~= j
                phi_jk = F/(R*T)*(Va(j) - Va(k));
                if abs(phi_jk) > 1e-10
                    IK_gap(j) = IK_gap(j) + G(P_K_gap,phi_jk,K_iA(j),K_iA(k));
                    INa_gap(j) = INa_gap(j) + G(P_Na_gap,phi_jk,Na_iA(j),Na_iA(k));
                else
                    IK_gap(j) = IK_gap(j) + P_K_gap*F*(K_iA(j) - K_iA(k));
                    INa_gap(j) = INa_gap(j) + P_Na_gap*F*(Na_iA(j) - Na_iA(k));
                end
            end
        end
    end
    Igap = IK_gap + INa_gap;
    
    % Dirichlet boundary conditions on extracellular potassium/sodium
    % around first/last neurons
    K_e0 = 3.5; K_eN = K_e0;
    Na_e0 = 135; Na_eN = Na_e0;
    
    % Injection of K+ to the extracellular space of the middle cells
    Iinj = zeros(N,1);
    mid = [24,25,26,27];
    if pumpK
        Iinj(mid) = eps;
    end
    
    % Model for a neuron-astrocyte network responsible for initiation and 
    % propagation of cortical spreading depolarization
    dx = zeros(10*N,1);
    dx(1:N) = -(1/Cm)*(INa + INaP + IK + IL + IPn);
    dx(N+1:2*N) = phi_n*(nInf - n)./theta_n;
    dx(2*N+1:3*N) = phi_hp*(hpInf - hp)./theta_hp;
    
    dx(3*N+1) = DNa*(Na_e(2) - 2*Na_e(1) + Na_e0)/delta_x^2 ...
                + 10*Sn/(F*omega_e)*(INa(1) + INaP(1) + 3*IPn(1)) ...
                + 10*Sa/(F*omega_e)*(INa_A(1) + 3*IPa(1));
    dx(3*N+2:4*N-1) = DNa*(Na_e(3:N) - 2*Na_e(2:N-1) + Na_e(1:N-2))/delta_x^2 ...
                      + 10*Sn/(F*omega_e)*(INa(2:N-1) + INaP(2:N-1) + 3*IPn(2:N-1)) ...
                      + 10*Sa/(F*omega_e)*(INa_A(2:N-1) + 3*IPa(2:N-1));
    dx(4*N) = DNa*(Na_eN - 2*Na_e(N) + Na_e(N-1))/delta_x^2 ...
              + 10*Sn/(F*omega_e)*(INa(N) + INaP(N) + 3*IPn(N)) ...
              + 10*Sa/(F*omega_e)*(INa_A(N) + 3*IPa(N));
    
    dx(4*N+1:5*N) = -10*Sn/(F*omega_n)*(INa + INaP + 3*IPn);
    
    dx(5*N+1) = DK*(K_e(2) - 2*K_e(1) + K_e0)/delta_x^2 ...
                + 10*Sn/(F*omega_e)*(IK(1) - 2*IPn(1)) ...
                + 10*Sa/(F*omega_e)*(IK_A(1) - 2*IPa(1));
    dx(5*N+2:6*N-1) = DK*(K_e(3:N) - 2*K_e(2:N-1) + K_e(1:N-2))/delta_x^2 ...
                      + 10*Sn/(F*omega_e)*(IK(2:N-1) - 2*IPn(2:N-1)) + Iinj(2:N-1) ...
                      + 10*Sa/(F*omega_e)*(IK_A(2:N-1) - 2*IPa(2:N-1));
    dx(6*N) = DK*(K_eN - 2*K_e(N) + K_e(N-1))/delta_x^2 ...
              + 10*Sn/(F*omega_e)*(IK(N) - 2*IPn(N)) ...
              + 10*Sa/(F*omega_e)*(IK_A(N) - 2*IPa(N));
    
    dx(6*N+1:7*N) = -10*Sn/(F*omega_n)*(IK - 2*IPn);
    
    dx(7*N+1:8*N) = -(1/Cma)*(INa_A + IK_A + IPa + Igap);
    
    dx(8*N+1:9*N) = -10*Sa/(F*omega_a)*(INa_A + 3*IPa + INa_gap); % sign of INa_gap
     
    dx(9*N+1:10*N) = -10*Sa/(F*omega_a)*(IK_A - 2*IPa + IK_gap); % sign of IK_gap
    
    if pumpK
        if all(Vn(mid) > -30)
            t
            pumpK = false;
        end
    end
end