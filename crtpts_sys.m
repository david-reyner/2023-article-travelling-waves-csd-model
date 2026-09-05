function [F, DF] = crtpts_sys(x0, rhoN, rhoA, hp, conc)

    % Initial point
    Vn = x0(1); Va = x0(2); K_e = x0(3);
    
    % Phi_Vn and its partial derivatives
    [Phi_Vn, dVPhi_Vn, dKPhi_Vn] = diff_PhiVn([Vn, K_e], rhoN, hp, conc);
    
    % Phi_Va and its partial derivatives
    [Phi_Va, dVPhi_Va, dKPhi_Va] = diff_PhiVa([Va, K_e], rhoA, conc);
    
    % Phi_Ke and its partial derivatives
    [Phi_Ke, dVnPhi_Ke, dVaPhi_Ke, dKePhi_Ke] = diff_PhiKe(x0, rhoN, rhoA, conc);
    
    % Target function and its jacobian matrix
    F = [Phi_Vn; Phi_Va; Phi_Ke];
    DF = [dVPhi_Vn     0      dKPhi_Vn;
             0      dVPhi_Va  dKPhi_Va;
          dVnPhi_Ke dVaPhi_Ke dKePhi_Ke];
end
