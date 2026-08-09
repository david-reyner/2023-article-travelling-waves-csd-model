
close all; clear all; clc;

format long;

global pumpK
pumpK = true;

N = 50; % number of neurons in the network model

% Initial conditions for the neuron-astrocyte network model
x0 = ones(3*N,1);
x0(1:N) = -70; % Vn
x0(N+1:2*N) = 3.5; % Ke
x0(2*N+1:3*N) = -75.1757; % Va

hp0 = 0.9751; % fixing gating variable hp
Na_e = 135; % Extracellular concentration of sodium
Na_i0 = 3.5; % Initial intracellular concentration of sodium
Na_iA = 3.5; % Intracellular concentration of sodium
K_e = 3.5; % Extracellular concentration of potassium
K_i0 = 135; % Initial intracellular concentration of potassium
K_iA = 135; % Intracellular concentration of potassium

% Fixed and total initial concentrations
conc = zeros(6,1);
conc(1) = Na_e; conc(2) = Na_i0; conc(3) = Na_iA;
conc(4) = K_e; conc(5) = K_i0; conc(6) = K_iA;

% Maximal strength of neurons and astrocytes pump currents (microA/cm2)
rhoN = 5; rhoA = 5;

% Integrating a reduced version of the neuron-astrocyte network model
options_ode = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
t0 = 0; tf = 25000; % initial/final time of the integration
[t, y] = ode15s(@(t, x) red_nan_model(t, x, N, rhoN, rhoA, hp0, conc), ...
                                      [t0 tf], x0, options_ode);
                                  
% Save data
name_file = ['csd_nan_rlee_red_model', '_rhon', num2str(rhoN), ...
             '_rhoA', num2str(rhoA), '_Vn.txt'];
fclose(fopen(name_file, 'w+')); % open/create file discarding existing content

file = fopen(name_file, 'a');
formatOutput = strjoin(repmat({'%16.15f'}, 1, N+1));
for i = 1:size(y,1)
    fprintf(file, [formatOutput, '\r\n'], [t(i) y(i,1:N)]);
end
fclose(file);

% Plotting time evolution of Vn, Ke and Va for all neuron-astrocyte pairs
nf = {'Vn', 'Ke', 'Va'}; % auxiliar names to save figures
ttl = ["Evolution of neurons' voltage", ...
       "Evolution of [K^+]_e", ...
       "Evolution of astrocytes' voltage"]; % titles (1st figure)
ttl1 = ["Neurons membrane potential", ...
        "Extracellular K^+ concentration", ...
        "Astrocytes membrane potential"]; % titles (2nd figure)
ylb = ["V_N (mV)", "[K^+]_e (mM)", "V_A (mV)"]; % ylabels (1st figure)
cxs = [-70 20; 5 150; -80 0]; % caxis (2nd figure)
for i = 1:3
    % Time evolution of Vn, Va and Ke
    figure; hold on; grid on; box on; set(gca, 'FontSize', 13);
    xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); 
    ylabel(ylb(i), 'FontSize', 15, 'FontWeight', 'bold');
    title(ttl(i), 'FontSize', 15, 'FontWeight', 'bold');
    h = plot(t/1000, y(:,N*(i-1)+1:i*N), 'Linewidth', 1.5);
    set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines

    % -- Colorbar
    c = colorbar; colormap(parula(N)); c.FontSize = 13;
    c.TickLabels = cellfun(@num2str, num2cell(0:10:N), 'UniformOutput', false); 
    set(c.Label, 'String', '# cell', 'Fontsize', 15, 'Fontweight', 'bold');

    % Save figure
    name_fig = ['csd_nan_rlee_red_model', '_rhon', num2str(rhoN), ...
                '_', nf{i}, '_time.eps'];
    print(gcf, '-depsc', '-tiff', name_fig);
    
    % Generation and propagation of a travelling front for Vn, Va and Ke
    figure; hold on; box on; set(gca, 'FontSize', 13);
    xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
    ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
    title(ttl1(i), 'Fontsize', 15, 'FontWeight', 'bold');
    surf(1:N, t/1000, y(:,N*(i-1)+1:i*N)); shading interp; caxis(cxs(i,:));
    set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
    axpos = get(gca, 'Position'); % get current position
    set(gca, 'Position', axpos + [-0.035 0 0 0]);
    cb = colorbar; % add colorbar to the Poincare section plot
    set(cb, 'Position', cb.Position + [0.1 0 0 0]);
    set(cb, 'Ticks', cxs(i,:)); % set first/last ticks and skip intermediate ones on colorbar
    set(cb.Label, 'String', ylb{i}, 'Fontsize', 14, 'Fontweight', 'bold', ...
        'VerticalAlignment', 'baseline'); % set label to the colorbar
    
    if i == 2
        set(cb.Ruler, 'Exponent', 1); % set exponent 1 to the colorbar
        set(cb.Ruler.SecondaryLabel, 'Position', [1 154 0]); % location of exponent
    end
    
    % Save figure
    name_fig = ['csd_nan_rlee_red_model', '_rhon', num2str(rhoN), ...
        '_', nf{i}, '.eps'];
    print(gcf, '-depsc', '-tiff', name_fig);
end
return
%%

format long;

N = 50;
tf = 25000;

% Read data
name_file = ['csd_nan_rlee_red_model', '_rhon', num2str(rhoN), ...
             '_rhoA', num2str(rhoA), '_N', num2str(N), ...
             '_tf', num2str(tf),'_Vn.txt'];
dat = importdata(name_file); t = dat(:,1); y = dat(:,2:end);

% -------------------- Estimation of the wave speed c ------------------- %
V_N = y(t/1000 >= 0.05,1:N);
t = t(t/1000 >= 0.05);

ind1 = 10; ind2 = 20; % neuron indexes to estimate c
Vn1 = griddedInterpolant(V_N(:,ind1),t/1000); % linear interpolation of Vn(:,ind1) solution
Vn2 = griddedInterpolant(V_N(:,ind2),t/1000); % linear interpolation of Vn(:,ind2) solution
delta_t = Vn1(-20) - Vn2(-20);% time difference between solutions at V = -20
c = abs(ind2 - ind1)/delta_t % estimation of wave speed (cells per second)
c = c*60*0.044 % estimation of wave speed (mm per minut)
% c = 0.321377329906812; % wave speed (mm^2 per minut)
% D = 0.1176; % diffusion coefficient (mm^2 per minut)

figure; hold on; box on; set(gca, 'FontSize', 13);
cmap = parula(N); set(gca, 'Colororder', cmap([ind1, ind2],:));
xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); 
ylabel('V_N (mV)', 'FontSize', 15, 'FontWeight', 'bold');
title('Estimation of wave speed', 'FontSize', 15, 'FontWeight', 'bold');
plot(t/1000, V_N(:,[ind1 ind2]), 'Linewidth', 1.5); % plotting Vn1 and Vn2 solutions
plot([Vn1(-20), Vn2(-20)], [-20 -20], 'ko', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'k');
text(1, 22, '$\displaystyle vel \approx \frac{\textnormal{idx}_1 - \textnormal{idx}_2}{\Delta t} \,$cells/s', 'Fontsize', 16, 'Interpreter', 'Latex');
p1 = [Vn1(-20) -20]; p2 = [Vn2(-20) -20]; % initial/final point
xlm = xlim; ylm = ylim; % current axes limits
pos = get(gca, 'Position'); % current position
x = pos(1) + pos(3)/diff(xlm)*([p1(1) p2(1)] - xlm(1)); % x-coordinates of the arrow
y = pos(2) + pos(4)/diff(ylm)*([p1(2) p2(2)] - ylm(1)); % y-coordinates of the arrow
annotation('doublearrow', x, y); % doublearrow indicating the computed distance
text((p1(1)+p2(1))/2, p1(2), '\Deltat', 'Fontsize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % text label
legend({['V_{', num2str(ind1), '}'], ['V_{', num2str(ind2), '}']}, ...
    'Fontsize', 14, 'Fontweight', 'bold', 'Location', 'SouthEast');

% Save figure
name_fig = ['wave_speed_estimation', num2str(N), '_tf', num2str(tf), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Compute c from all possible pairs of neuron indices
ws = zeros(N, N);
for i = 1:N
    for j = 1:N
        if all(i+j ~= [2*i, N+1]) % remove diagonal combinations
            ind1 = i; ind2 = j; % neuron indexes to estimate c
            Vn1 = griddedInterpolant(V_N(:,ind1),t/1000); % linear interpolation of Vn(:,ind1) solution
            Vn2 = griddedInterpolant(V_N(:,ind2),t/1000); % linear interpolation of Vn(:,ind2) solution
            delta_t = Vn1(-20) - Vn2(-20);% time difference between solutions at V = -20
            ws(i,j) = (ind1 - ind2)/delta_t; % estimation of wave speed (cells per second)
            ws(i,j) = ws(i,j)*60*0.044; % estimation of wave speed (mm per minut)
        end
    end
end
ws(abs(ws) > 50) = nan; % remove wave speeds greater than -50
ws(abs(ws) == 0) = nan; % remove zero wave speeds

% Save data
name_file = ['wave_speeds_pairwise_rhon', num2str(rhoN), ...
             '_rhoA', num2str(rhoA), '_N', num2str(N), ...
             '_tf', num2str(tf), '.txt'];
fclose(fopen(name_file, 'w+'));

file = fopen(name_file, 'a');
formatOutput = strjoin(repmat({'%16.15f'}, 1, N));
for i = 1:size(ws,1)
    fprintf(file, [formatOutput, '\r\n'], ws(i,:));
end
fclose(file);

% Save figure
% figure; imagesc(ws); axis image;
figure; axis image; pc = pcolor(ws); pc.LineStyle = 'none';
xlabel('# cell', 'Fontsize', 15, 'Fontweight', 'bold');
ylabel('# cell', 'Fontsize', 15, 'Fontweight', 'bold');
set(gca, 'Fontsize', 13, 'Layer', 'top'); cb = colorbar; caxis([-50 50]);
set(cb.Label, 'String', 'c', 'Fontsize', 15, 'Fontweight', 'bold', ...
    'Rotation', 0);
title('Pairwise wave speeds', 'Fontsize', 15, 'Fontweight', 'bold');

name_fig = ['pairwise_wave_speeds_N', num2str(N), '_tf', num2str(tf), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

figure; hold on; box on; set(gca, 'Fontsize', 13);
xlabel('c', 'Fontsize', 15, 'Fontweight', 'bold');
ylabel('# cells', 'Fontsize', 15, 'Fontweight', 'bold');
hist(ws(:), -50:1:50);

name_fig = ['hist_wave_speeds_N', num2str(N), '_tf', num2str(tf), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);
saveas(gcf, ['hist_wave_speeds_N', num2str(N), '_tf', num2str(tf), '.fig']);
%-------------------------------------------------------------------------%