
close all; clear all; clc;

format long;

global pumpK
pumpK = true;

N = 50; % number of neurons in the network model
delta_x = 1; % distance between neurons

% Initial conditions for the neuron-astrocyte network model
x0 = ones(10*N,1);
x0(1:N) = -70; % Vn
x0(N+1:2*N) = 0.25512; % n
x0(2*N+1:3*N) = 0.9751; % hp
x0(3*N+1:4*N) = 133.574; % Na_e
x0(4*N+1:5*N) = 4.297; % Na_i
x0(5*N+1:6*N) = 3.5; % K_e
x0(6*N+1:7*N) = 138.116; % K_i
x0(7*N+1:8*N) = -75.1757; % Va
x0(8*N+1:9*N) = 6.6; % Na_iA
x0(9*N+1:10*N) = 124; % K_iA

% Rate of injection of K+ into the extracellular space
eps = 0.005;

% Maximal strength of neurons and astrocytes pump currents (microA/cm2)
rhoN = 5; rhoA = 5;

% Maximal conductance of a gap junction
sig_gap = 0.0;

% Number of gap junction connections
Ncn = 3;

% Integrating neuron-astrocyte network model
options_ode = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
t0 = 0; tf = 70000; % initial/final time of the integration
[t, y] = ode15s(@(t, x) neuron_astrocyte_RLee_model(t, x, delta_x, N, ...
                        eps, rhoN, rhoA, sig_gap, Ncn), [t0 tf], x0);

% V(x,t) vs t
figure; hold on; grid on; box on; set(gca, 'FontSize', 13);
xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); xlim([t0 tf]/1000);
ylabel('V (mV)', 'FontSize', 15, 'FontWeight', 'bold');
title('Time evolution of neurons voltage', 'FontSize', 15, 'FontWeight', 'bold');
h = plot(t/1000, y(:,1:N), 'Linewidth', 1.5);
set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines

% -- Colorbar
c = colorbar; c.FontSize = 13;
c.TickLabels = cellfun(@num2str, num2cell(0:10:N), 'UniformOutput', false); 
set(c.Label, 'String', '# cell', 'Fontsize', 15, 'Fontweight', 'bold');

% Save figure
name_fig = ['csd_nan_rlee_VnTime_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_sgap', num2str(sig_gap), '_N', num2str(Ncn), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% V(x,t) vs x
% figure; hold on; box on; set(gca, 'FontSize', 13);
% xlabel('x', 'FontSize', 15, 'FontWeight', 'bold'); 
% ylabel('V (mV)', 'FontSize', 15, 'FontWeight', 'bold');
% title('Spatial evolution of neurons voltage', 'FontSize', 15, 'FontWeight', 'bold');
% h = plot(1:N, y(:,1:N), 'Linewidth', 1.5);
% set(h, {'Color'}, num2cell(parula(length(t)),2)); % change default color of lines

% Neurons' membrane potential
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Neurons membrane potential', 'Fontsize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, y(:,1:N)); shading interp; caxis([-70 -10]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', 'V_N (mV)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar

% Save figure
name_fig = ['csd_nan_rlee_Vn_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_sgap', num2str(sig_gap), '_N', num2str(Ncn), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Astrocytes' membrane potential
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Astrocytes membrane potential', 'Fontsize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, y(:,7*N+1:8*N)); shading interp; caxis([-80 -20]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', 'V_A (mV)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar

% Save figure
name_fig = ['csd_nan_rlee_Va_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_sgap', num2str(sig_gap), '_N', num2str(Ncn), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Extracellular potassium concentration
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Extracellular K^+ concentration', 'Fontsize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, y(:,5*N+1:6*N)); shading interp; caxis([3.5 70]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', '[K^+]_e (mM)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar

% Save figure
name_fig = ['csd_nan_rlee_K_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_sgap', num2str(sig_gap), '_N', num2str(Ncn), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Extracellular sodium concentration
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Extracellular Na^+ concentration', 'Fontsize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, y(:,3*N+1:4*N)); shading interp; caxis([20 140]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', '[Na^+]_e (mM)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'bottom'); % set label to the colorbar

% Save figure
name_fig = ['csd_nan_rlee_Na_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_sgap', num2str(sig_gap), '_N', num2str(Ncn), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Plotting time evolution of the variables for a fixed neuron/astrocyte pair
ind = 15; % neuron/astrocyte index
k = [0 7; 3 5; 4 6]; % pairs of the indexes of the variables to plot
lg = {'V_N', 'V_A'; '[Na^+]_e', '[K^+]_e'; '[Na^+]_i', '[K^+]_i'}; % legend labels
pos = ["NorthWest", "West", "West"]; % legend locations
ylb = ["Neuron/astrocyte membrane potential", ...
       "Extracellular Na^+, K^+ concentrations", ...
       "Intracellular Na^+, K^+ concentrations"]; % ylabels
for i = 1:3
    figure; hold on; box on; set(gca, 'FontSize', 13);
    xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel(ylb(i), 'FontSize', 15, 'FontWeight', 'bold');
    plot(t/1000, y(:,k(i,1)*N+ind), 'b-', 'Linewidth', 1.5);
    plot(t/1000, y(:,k(i,2)*N+ind), 'r-', 'Linewidth', 1.5);
    legend(lg{i,:}, 'Fontsize', 15, 'Fontweight', 'bold', 'Location', pos(i));
end