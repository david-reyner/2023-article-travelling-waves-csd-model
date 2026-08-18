
close all; clear all; clc;

format long;

% Initial sodium/potassium concentrations
hp0 = 0.9751; % fixing gating variable hp
Na_e = 135; % Extracellular concentration of sodium
Na_i0 = 3.5; % Initial intracellular concentration of sodium
Na_iA = 3.5; % Intracellular concentration of sodium
K_e0 = 3.5; % Extracellular concentration of potassium
K_i0 = 135; % Initial intracellular concentration of potassium
K_iA = 135; % Intracellular concentration of potassium

% Fixed and total initial concentrations
conc = zeros(6,1);
conc(1) = Na_e; conc(2) = Na_i0; conc(3) = Na_iA;
conc(4) = K_e0; conc(5) = K_i0; conc(6) = K_iA;

% Maximal strength of neurons and astrocytes pump currents (microA/cm2)
rhoN = 5; rhoA = 5;

%%%%%%%%%%%%% 3D plots of functions Phi_Vn, Phi_Va and Phi_Ke %%%%%%%%%%%%%
K_e = linspace(1,160); V_n = linspace(-100,20); V_a = linspace(-100,20);

% Plotting Phi_Vn as a function of V_n and K_e
Phi_Vn = zeros(length(V_n), length(K_e));
for i = 1:length(V_n)
    for j = 1:length(K_e)
        Phi_Vn(i,j) = rhs_red_nan(V_n(i), 0, K_e(j), rhoN, rhoA, hp0, conc, 1, 1);
    end
end

figure; hold on; grid on; box on; set(gca, 'Fontsize', 13);
xlabel('$\mathbf{z}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex'); xlim([K_e(1) K_e(end)]);
ylabel('$\mathbf{x}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex'); ylim([V_n(1) V_n(end)]);
zlabel('$\mathbf{\widetilde{f}(x,z)}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex');
surf(K_e, V_n, Phi_Vn); shading flat;
surf(K_e, V_n, zeros(size(Phi_Vn)), 'FaceColor', [0.3 0.3 0.3], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
set(gca, 'Position', get(gca, 'Position') + [0.035 0 0 0]);
view([-40 30]);

name_fig = ['fxz_rhon', num2str(rhoN), '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Plotting Phi_Va as a function of V_a and K_e
Phi_Va = zeros(length(V_a), length(K_e));
for i = 1:length(V_a)
    for j = 1:length(K_e)
        Phi_Va(i,j) = rhs_red_nan(0, V_a(i), K_e(j), rhoN, rhoA, hp0, conc, 2, 1);
    end
end

figure; hold on; grid on; box on; set(gca, 'Fontsize', 13);
xlabel('$\mathbf{z}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex'); xlim([K_e(1) K_e(end)]);
ylabel('$\mathbf{y}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex'); ylim([V_a(1) V_a(end)]);
zlabel('$\mathbf{\widetilde{g}(y,z)}$', 'FontSize', 16, 'FontWeight', 'bold', ...
       'Interpreter', 'Latex');
surf(K_e, V_a, Phi_Va); shading flat;
surf(K_e, V_a, zeros(size(Phi_Va)), 'FaceColor', [0.3 0.3 0.3], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
view([-40 30]);

name_fig = ['gyz_rhon', num2str(rhoN), '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Plotting Phi_Ke as a function of V_n and K_e (for different values V_a)
y = [-80 -20 0 50];
Phi_Ke = zeros(length(V_n), length(K_e));
for k = 1:length(y)
    for i = 1:length(V_n)
        for j = 1:length(K_e)
            Phi_Ke(i,j) = rhs_red_nan(V_n(i), y(k), K_e(j), rhoN, rhoA, hp0, conc, 3, 1);
        end
    end

    figure; hold on; grid on; box on; set(gca, 'Fontsize', 13);
    xlabel('$\mathbf{z}$', 'FontSize', 16, 'FontWeight', 'bold', ...
           'Interpreter', 'Latex'); xlim([K_e(1) K_e(end)]);
    ylabel('$\mathbf{x}$', 'FontSize', 16, 'FontWeight', 'bold', ...
           'Interpreter', 'Latex'); ylim([V_n(1) V_n(end)]);
    zlabel('$\mathbf{h(x,y,z)}$', 'FontSize', 15, 'FontWeight', 'bold', ...
           'Interpreter', 'Latex');
    title(['$\mathbf{y = ', num2str(y(k)), '}$'], 'FontSize', 15, 'FontWeight', 'bold', ...
           'Interpreter', 'Latex');
    surf(K_e, V_n, Phi_Ke); shading flat;
    surf(K_e, V_n, zeros(size(Phi_Vn)), 'FaceColor', [0.3 0.3 0.3], ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none');
    view([-40 30]);

    name_fig = ['hxyz_y', num2str(y(k)), '_rhon', num2str(rhoN), '_rhoa', num2str(rhoA), '.eps'];
    print(gcf, '-depsc', '-tiff', name_fig);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%