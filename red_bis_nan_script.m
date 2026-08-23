
close all; clear all; clc;

format long;

global pumpK
pumpK = true;

% Setting parameter values
N = 50; % number of neurons in the network model

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

% Storing parameters
parm = [rhoN, rhoA, hp0, conc'];

% Initial conditions to integrate a reduced neuron-astrocyte network model
x0 = ones(3*N,1);
x0(1:N) = -70; % Vn
x0(N+1:2*N) = 3.5; % Ke
x0(2*N+1:3*N) = -75.1757; % Va

% Integrating a reduced version of the neuron-astrocyte network model (not
% constrained to the critical manifold)
options_ode = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
t0 = 0; tf = 25000; % initial/final time of the integration
[t, y] = ode15s(@(t, x) red_nan_model(t, x, N, rhoN, rhoA, hp0, conc), ...
                                      [t0 tf], x0, options_ode);

%%%%%%%%%%% Computation of the zero curves of Phi_Vn and Phi_Va %%%%%%%%%%%
Vn0 = -70; % initial seed for the zero curve of Phi_Vn
Va0 = -75; % initial seed for the zero curve of Phi_Va
Ke0 = 1; % initial parameter value

% Continuate the zero curve of Phi_Va (no folds)
Va0c = cont_PhiVa(Va0, rhoA, conc);

% Compute an initial solution of the zero curve of Phi_Vn
[Vn0, ~] = alg_constraints(Vn0, Va0, Ke0, rhoN, rhoA, hp0, conc);
x0 = [Vn0; Ke0];

% Continuate initial solution to find the zero curve of Phi_Vn
Vn0c = cont_PhiVn(x0, 1, rhoN, hp0, conc);
Vn0cbk = cont_PhiVn(x0, -1, rhoN, hp0, conc);
Vn0c = [flipud(Vn0cbk(2:end,:)); Vn0c];

% Continuate initial solution to find the zero curve of Phi_Vn: Ke(Vn)
Ke0c = cont_PhiVn2([Ke0; Vn0], 1, rhoN, hp0, conc);
Ke0cbk = cont_PhiVn2([Ke0; Vn0], -1, rhoN, hp0, conc);
Ke0c = [flipud(Ke0cbk(2:end,:)); Ke0c];

% Linear interpolation of zero curves Ke(Vn) and Va(Ke)
Ke0cp = griddedInterpolant(Ke0c(:,2), Ke0c(:,1), 'linear', 'none'); % Ke(Vn) curve
Va0cp = griddedInterpolant(Va0c(:,2), Va0c(:,1), 'linear', 'none'); % Va(Ke) curve

% Zero curve of Phi_Vn: Ke(Vn)
figure; hold on; box on; set(gca, 'Fontsize', 13);
xlabel('V_N', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('[K^+]_e', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 160]);
title('\Phi_{V_N} = 0', 'FontSize', 13, 'FontWeight', 'bold');
plot(Ke0c(:,2), Ke0c(:,1), '-', 'Color', [0.8 0.15 0.15], 'Linewidth', 1.5);

% Zero curves of Phi_Vn and Phi_Va
fig_zr_crvs = figure; hold on; box on; set(gca, 'Fontsize', 13); ax_pr = gca;
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([0 160]);
ylabel('X^\ast(z) , Y(z)', 'FontSize', 15, 'FontWeight', 'bold');
title('f(x,z) = g(y,z) = 0', 'FontSize', 15, 'FontWeight', 'bold');
plot(Vn0c(:,2), Vn0c(:,1), '-', 'Color', [0.8 0.15 0.15], 'Linewidth', 2);
plot(Va0c(:,2), Va0c(:,1), '-', 'Color', [1 0.5 0], 'Linewidth', 2);
legend({'x', 'y'}, 'Fontsize', 15, 'Fontweight', 'bold', ...
    'Location', 'NorthWest', 'Box', 'off', 'AutoUpdate', 'off');
text(ax_pr, 80, 10, '$X^r$', 'Fontsize', 18, 'Interpreter', 'latex');

% -- Create an inset figure
xstart = 0.55; xend = 0.85; ystart = 0.23; yend = 0.63;
ax_ch = axes('Position',[xstart, ystart, xend - xstart, yend - ystart]);
set(gca, 'Fontsize', 13, 'XLim', [0 30], 'yLim', [-90 -10]); hold on; box on;

% ------------------------------ Save data ------------------------------ %
% --> Phi_Vn
name_file = ['zero_curves_phiVn_rhon', num2str(rhoN), ...
             '_rhoa', num2str(rhoA), '.txt'];
fclose(fopen(name_file, 'w+'));

file = fopen(name_file, 'a');
for i = 1:size(Vn0c,1)
    fprintf(file, '%16.15f %16.15f %16.15f %d\r\n', Vn0c(i,:));
end
fclose(file);

% --> Phi_Va
name_file = ['zero_curves_phiVa_rhon', num2str(rhoN), ...
             '_rhoa', num2str(rhoA), '.txt'];
fclose(fopen(name_file, 'w+'));

file = fopen(name_file, 'a');
for i = 1:size(Va0c,1) 
    fprintf(file, '%16.15f %16.15f %16.15f %d\r\n', Va0c(i,:));
end
fclose(file);
%-------------------------------------------------------------------------%

% Determining first/second folds of the curve Vn0c (zeros of dPhi_Vn/dVn)
% Ke = @(x) PhiVn0(Ke0cp(x), x, rhoN, hp0, conc);
df1 = @(x) diff_PhiVn([x, Ke0cp(x)], rhoN, hp0, conc, 2); % dPhi_Vn/dVn
Vn = df1(Ke0c(:,2));
ind1 = find(Vn(1:end-1) < 0 & Vn(2:end) > 0);
ind2 = find(Vn(1:end-1) > 0 & Vn(2:end) < 0);
fld_Vn = secant_method(df1, Ke0c(ind1,2), Ke0c(ind1+1,2), 5e-12); % 1st fold
fld_Ke = Ke0cp(fld_Vn);
fld_Vn2 = secant_method(df1, Ke0c(ind2,2), Ke0c(ind2+1,2), 5e-12); % 2nd fold
fld_Ke2 = Ke0cp(fld_Vn2);

ind1 = find(Vn0c(:,1) < fld_Vn); % up to 1st fold
ind2 = find(Vn0c(:,2) >= fld_Ke); % upper branch starting after 1st fold

indK1 = find(Ke0c(:,2) <= fld_Vn); % up to maximum
indK2 = find(Ke0c(:,1) > fld_Ke); % rightmost branch above maximum

% Indexes for each branch of Vn0c
indL = find(Vn0c(:,1) <= fld_Vn); % lower branch's indexes
indM = find(Vn0c(:,1) > fld_Vn & Vn0c(:,1) < fld_Vn2); % middle branch's indexes
indM = [indM(1)-1; indM; indM(end)+1]; % include indices before/after folds
indM = flipud(indM); % auxiliar step for subsequent interpolation with griddedInterpolant
indU = find(Vn0c(:,1) >= fld_Vn2); % upper branch's indexes

% Linear interpolation of zero curves Phi_Vn (3 branches) and Phi_Va
Vn0cpL = griddedInterpolant(Vn0c(indL,2), Vn0c(indL,1), 'linear', 'none'); % lower branch
Vn0cpM = griddedInterpolant(Vn0c(indM,2), Vn0c(indM,1), 'linear', 'none'); % middle branch
Vn0cpU = griddedInterpolant(Vn0c(indU,2), Vn0c(indU,1), 'linear', 'none'); % upper branch
Vn0cp = griddedInterpolant(Vn0c(unique([ind1; ind2]),2), ...
            Vn0c(unique([ind1; ind2]),1), 'linear', 'none'); % Vn(Ke) curve

% -- Plot folds in the curve Vn0c
plot(ax_pr, [fld_Ke, fld_Ke], [-120, fld_Vn], ':', 'Linewidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
plot(ax_pr, [fld_Ke2, fld_Ke2], [-120, fld_Vn2], ':', 'Linewidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
plot(ax_pr, [0, fld_Ke], [fld_Vn, fld_Vn], ':', 'Linewidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
plot(ax_pr, [0, fld_Ke2], [fld_Vn2, fld_Vn2], ':', 'Linewidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
plot(ax_pr, fld_Ke, fld_Vn, 'ko', 'Linewidth', 2, 'MarkerFaceColor', 'k');
plot(ax_pr, fld_Ke2, fld_Vn2, 'ko', 'Linewidth', 2, 'MarkerFaceColor', 'k');

% -- Copy all child objects to inset figure
copyobj(allchild(ax_pr), ax_ch);
set(ax_ch, 'XTick', [0 fld_Ke2 fld_Ke, 30], 'XTickLabel', {'0', 'z^L', 'z^R', '30'})
set(ax_ch, 'YTick', [-80 fld_Vn fld_Vn2, -10], 'YTickLabel', {'-80', 'x^R', 'x^L', '-10'})
text(ax_ch, (fld_Ke2+fld_Ke)/2-2, -75, '$X^l$', 'Fontsize', 18, 'Interpreter', 'latex');
text(ax_ch, fld_Ke-3.25, (fld_Vn2+fld_Vn)/2+4, '$X^m$', 'Fontsize', 18, 'Interpreter', 'latex');

% -- Draw a rectangle in parent figure to indicate which part has been zoomed in
rectangle(ax_pr, 'Position', [0 -80 30 70], 'EdgeColor', 'k', ...
          'Linestyle', '-', 'Linewidth', 0.5);

% -- Join rectangle's left-bottom point with inset's left-bottom one
pos_pr = ax_pr.Position; pos_ch = ax_ch.Position;
x_pr = pos_pr(1) + pos_pr(3)*(0 - ax_pr.XLim(1))/diff(ax_pr.XLim);
y_pr = pos_pr(2) + pos_pr(4)*(-80 - ax_pr.YLim(1))/diff(ax_pr.YLim);
x_ch = pos_ch(1);
y_ch = pos_ch(2);
annotation(fig_zr_crvs, 'line', [x_pr x_ch], [y_pr y_ch], 'Color', 'k', ...
           'Linestyle', '--', 'Linewidth', 1)

% -- Join rectangle's right-top point with inset's right-top one
pos_pr = ax_pr.Position; pos_ch = ax_ch.Position;
x_pr = pos_pr(1) + pos_pr(3)*(30 - ax_pr.XLim(1))/diff(ax_pr.XLim);
y_pr = pos_pr(2) + pos_pr(4)*(-10 - ax_pr.YLim(1))/diff(ax_pr.YLim);
x_ch = pos_ch(1) + pos_ch(3);
y_ch = pos_ch(2) + pos_ch(4);
annotation(fig_zr_crvs, 'line', [x_pr x_ch], [y_pr y_ch], 'Color', 'k', ...
           'Linestyle', '--', 'Linewidth', 1)

% Save figure
name_fig = ['zero_curves_phiVn_phiVa_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);
       
% 3D curve defined by the two algebraic constraints (critical manifold)
figure; hold on; box on; grid on; set(gca, 'Fontsize', 13);
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([0 160]);
ylabel('x', 'FontSize', 15, 'FontWeight', 'bold');
zlabel('y', 'FontSize', 15, 'FontWeight', 'bold');
title('\textbf{Critical manifold} $\mathbf{\mathcal{\widetilde{M}}_0}$', ...
      'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'Latex');
plot3(Vn0c(:,2), Vn0c(:,1), Va0cp(Vn0c(:,2)), '--', 'Color', [0.3 0.3 0.3], ...
    'Linewidth', 2);
% legend('$\mathcal{\widetilde{M}}_0$', 'Fontsize', 14, 'Box', 'off', 'Fontweight', 'bold', ...
%        'Location', 'West', 'Interpreter', 'Latex');
view([-30,30]);

% Plot solution of the reduced model (for a fixed neuron-astrocyte pair)
indna = 20; cmap = parula(N);
plot3(y(:,N+indna), y(:,indna), y(:,2*N+indna), '-', ...
    'Linewidth', 2, 'Color', cmap(indna,:));
% set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines
plot3(y(1,N+indna), y(1,indna), y(1,2*N+indna), 'ko', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k', 'Linewidth', 1.5); % first point

% Save figure
name_fig = ['critical_manifold_phiVn_phiVa_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Critical manifold
figure; hold on; box on; grid on; set(gca, 'Fontsize', 13);
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 275]);
ylabel('w', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 1.5]);
zlabel('x, y', 'FontSize', 15, 'FontWeight', 'bold'); zlim([-80 50]);
Z = repmat(Vn0c(:,2), [1 2]); Z1 = repmat(Va0c(:,2), [1 2]);
W = [0 1.5].*ones(length(Vn0c),1); W1 = [0 1.5].*ones(length(Va0c),1);
X = repmat(Vn0c(:,1), [1 2]); Y = repmat(Va0c(:,1), [1 2]);
surf(Z, W, X, 'FaceColor', [0.8 0.15 0.15], 'FaceAlpha', 0.5, ...
    'MeshStyle', 'col', 'EdgeColor', [0.8 0.15 0.15], ...
    'Linewidth', 1.5); % zero curve Phi_Vn time w axis
surf(Z1, W1, Y, 'FaceColor', [1 0.5 0], 'FaceAlpha', 0.5, ...
    'MeshStyle', 'col', 'EdgeColor', [1 0.5 0], ...
    'Linewidth', 1.5); % zero curve Phi_Va time w axis
view([-30,30]);

% Save figure
name_fig = ['zero_surfaces_phiVn_phiVa_rhon', num2str(rhoN), ...
            '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% --> Critical manifold: zero curve Phi_Va time w axis
% figure; hold on; box on; grid on; set(gca, 'Fontsize', 13);
% xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 275]);
% ylabel('w', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 1.5]);
% zlabel('y', 'FontSize', 15, 'FontWeight', 'bold'); zlim([-80 50]);
% Z = repmat(Va0c(:,2), [1 2]);
% W = [0 1.5].*ones(length(Va0c),1);
% Y = repmat(Va0c(:,1), [1 2]);
% surf(Z, W, Y, 'FaceColor', [1 0.5 0], 'FaceAlpha', 0.5, ...
%     'MeshStyle', 'col', 'EdgeColor', [1 0.5 0], 'Linewidth', 1.5);
% view([-30,30]);

% ------------------------------ Save data ------------------------------ %
% Saving data of the critical manifold
name_file = 'critical_manifold.txt';
fclose(fopen(name_file, 'w'));
file = fopen(name_file, 'a');
for i = 1:length(Vn0c(:,1))
    fprintf(file, '%16.15f %16.15f %16.15f\r\n', Vn0c(i,1), Vn0c(i,2), ...
        Va0cp(Vn0c(i,2)));
end
fclose(file);

% Saving data of the zero curves
name_file = 'zero_curve_Vn.txt';
fclose(fopen(name_file, 'w'));
file = fopen(name_file, 'a');
for i = 1:length(Vn0c(:,1))
    fprintf(file, '%16.15f %16.15f\r\n', Vn0c(i,1), Vn0c(i,2));
end
fclose(file);

name_file = 'zero_curve_Va.txt';
fclose(fopen(name_file, 'w'));
file = fopen(name_file, 'a');
for i = 1:length(Va0c(:,1))
    fprintf(file, '%16.15f %16.15f\r\n', Va0c(i,1), Va0c(i,2));
end
fclose(file);
%-------------------------------------------------------------------------%

% ------------------ Eigenvalues of the fast subsystem ------------------ %
% Stability of the critical manifold branches
% Fast subsystem: d/dt V_N = (1/c) \Phi_{V_N}(V_N, [K^+]_e)
%                 d/dt V_A = (1/c) \Phi_{V_A}(V_A, [K^+]_e)
c = 1; % 7.3 % -2.766678115588945; % value of the wave speed (cell/s or mm/min)
[~, lam1, ~] = diff_PhiVn(Vn0c(:,1:2), rhoN, hp0, conc);
[~, lam2, ~] = diff_PhiVa(Va0c(:,1:2), rhoA, conc);

% Set colors for each branch
colors = {[0 0 1], [1 0 1], [1 0 0]};
% colors = {[0, 0.447, 0.741], [0.494, 0.184, 0.556], [0.6350 0.0780 0.1840]};

% Plotting eigenvalues of the fast subsystem
xlims = {[1, 160], [1 25]};

% Eigenvalues of the fast susbsytem on different axes
figure; hold on; box on; grid on; set(gca, 'Fontsize', 13);
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim(xlims{1});
title("Eigenvalues of the fast subsystem", 'Fontsize', 14, 'Fontweight', 'bold');

yyaxis right;
ylabel('\lambda_2', 'FontSize', 15, 'FontWeight', 'bold');
p1 = plot(Va0c(:,2), bsxfun(@times, (1./c), lam2), '-', 'Linewidth', 2); % lambda 2

yyaxis left; set(gca, 'YColor', 'k');
ylabel('\lambda_1^\ast', 'FontSize', 15, 'FontWeight', 'bold');
xline(fld_Ke, 'k--', 'Linewidth', 1.5); % 1st fold
xline(fld_Ke2, 'k--', 'Linewidth', 1.5); % 2nd fold
p2 = plot(Vn0c(indL,2), bsxfun(@times, (1./c), lam1(indL)), '-', ...
    'Color', colors{1}, 'Linewidth', 2); % lambda 1 lower branch
plot(Vn0c(indM,2), bsxfun(@times, (1./c), lam1(indM)), '-', ...
    'Color', colors{2}, 'Linewidth', 2); % lambda 1 middle branch
plot(Vn0c(indU,2), bsxfun(@times, (1./c), lam1(indU)), '-', ...
    'Color', colors{3}, 'Linewidth', 2); % lambda 1 right branch

% Create an inset figure
xstart = 0.42; xend = 0.77; ystart = 0.515; yend = 0.89;
axes('Position',[xstart ystart xend - xstart yend - ystart]);

hold on; box on; grid on; set(gca, 'Fontsize', 13); xlim(xlims{2});
yyaxis right;
plot(Va0c(:,2), bsxfun(@times, (1./c), lam2), '-', 'Linewidth', 2); % lambda 2

yyaxis left; set(gca, 'YColor', 'k'); set(gca, 'YLim', [-40 25]);
xline(fld_Ke, 'k--', 'Linewidth', 1.5); % 1st fold
xline(fld_Ke2, 'k--', 'Linewidth', 1.5); % 2nd fold
plot(Vn0c(indL,2), bsxfun(@times, (1./c), lam1(indL)), '-', ...
    'Color', colors{1}, 'Linewidth', 2); % lambda 1 lower branch
plot(Vn0c(indM,2), bsxfun(@times, (1./c), lam1(indM)), '-', ...
    'Color', colors{2}, 'Linewidth', 2); % lambda 1 middle branch
plot(Vn0c(indU,2), bsxfun(@times, (1./c), lam1(indU)), '-', ...
    'Color', colors{3}, 'Linewidth', 2); % lambda 1 right branch

% Line transparency
trns = [vertcat(p1.Color) 1]; % cs lists to matrix
col_trns = mat2cell(trns, ones(1,length(c))); % conversion matrix to cell
[p1.Color] = col_trns{:,1}; % assign variables by cs lists

trns = [vertcat(p2.Color) 1]; % cs lists to matrix [0.4 0.55 0.7 0.85 1]
col_trns = mat2cell(trns, ones(1,length(c))); % conversion matrix to cell
[p2.Color] = col_trns{:,1}; % assign variables by cs lists

% Save figure
name_fig = ['stability_eqpoints_fast_subsystem_rhon', num2str(rhoN), ...
            '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%% Equilibrium points of the reduced neuron-astrocyte network model %%%%
%%%%%%%%%%%%%%%% (equilibrium points of the slow subsystem) %%%%%%%%%%%%%%%

% Plotting zero level surfaces of Phi_Vn, Phi_Va and Phi_Ke
figure; hold on; box on; grid on; set(gca, 'Fontsize', 13); view([-37.5 30]);
xlabel('V_N', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('V_A', 'FontSize', 15, 'FontWeight', 'bold');
zlabel('[K^+]_e', 'FontSize', 15, 'FontWeight', 'bold');
title('\Phi_{V_N} = 0 , \Phi_{V_A} = 0 , \Phi_{[K^+]_e} = 0', ...
    'FontSize', 13, 'FontWeight', 'bold');

PhiVn = @(x,y,z) rhs_red_nan(x, y, z, rhoN, rhoA, hp0, conc, 1, 1); % get Phi_Vn
PhiVa = @(x,y,z) rhs_red_nan(x, y, z, rhoN, rhoA, hp0, conc, 2, 1); % get Phi_Va
PhiKe = @(x,y,z) rhs_red_nan(x, y, z, rhoN, rhoA, hp0, conc, 3, 1); % get Phi_Ke
int = [-80 80 -100 50 1 400]; % interval where to plot zero levels for each function
fimplicit3(PhiVn, int, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.5); % plot Phi_Vn = 0
fimplicit3(gca, PhiVa, int, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5); % plot Phi_Va = 0
fimplicit3(gca, PhiKe, int, 'EdgeColor', 'none', 'FaceColor', 'g'); % plot Phi_Ke = 0
axf3 = gca;

% Add critical manifold to the fimplicit plot
Ke = linspace(1,400,4000);
plot3(axf3, Vn0cp(Ke), Va0cp(Ke), Ke, 'k-', 'Linewidth', 1.5);

% ----------- Function Phi_Ke restricted to critical manifold ----------- %
% Plotting Phi_Ke restricted to the zero curves of Phi_Vn and Phi_Va

% Case 1: consider all the three branches of Phi_Vn
Kend = 300;
np = 3250;
Ke = [linspace(0.0001,fld_Ke,np), ... % discretization for the lower branch of Phi_Vn
      fliplr(linspace(fld_Ke2,fld_Ke,np)), ... % discretization for the middle branch of Phi_Vn
      linspace(fld_Ke2,Kend,np)]'; % discretization for the upper branch of Phi_Vn
Phi_Ke = zeros(length(Ke),1);

% Compute function Phi_Ke on the whole curve Phi_Vn = 0 (on the 3 branches)
for i = 1:length(Ke)
    if i <= np
        Vn0cpBr = Vn0cpL; % take the lower branch of Vn0cp
    elseif i <= 2*np
        Vn0cpBr = Vn0cpM; % take the middle branch of Vn0cp
    else
        Vn0cpBr = Vn0cpU; % take the upper branch of Vn0cp
    end
    Phi_Ke(i) = rhs_red_nan(Vn0cpBr(Ke(i)), Va0cp(Ke(i)), Ke(i), rhoN, ...
                            rhoA, hp0, conc, 3);
end

% Remove duplicates
% Ke([np+1,2*np]) = []; Phi_Ke([np+1,2*np]) = [];

% Set values at folds
Phi_Ke([np+1,2*np]) = Phi_Ke([np,2*np+1]);

% Plot Phi_Ke restricted to the whole critical manifold
figure; hold on; box on; set(gca, 'Fontsize', 13);
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([Ke(1) 250]);
ylabel('H^\ast(z)', 'FontSize', 15, 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'Exponent', -2); % set exponent to -2
plot(Ke(1:np), Phi_Ke(1:np), '-', ...
    'Color', colors{1}, 'Linewidth', 2); % lower branch
plot(Ke(np:2*np+1), Phi_Ke(np:2*np+1), '-', ...
    'Color', colors{2}, 'Linewidth', 2); % middle branch
plot(Ke(2*np+1:end), Phi_Ke(2*np+1:end), '-', ...
    'Color', colors{3}, 'Linewidth', 2); % right branch
plot([Ke(1) Ke(end)], [0 0], 'k--', 'Linewidth', 2); % zero line
Phi_Ke1 = Phi_Ke; Ke1 = Ke;

% Case 2: consider the lower branch and a piece of the third one of Phi_Vn
Ke = [linspace(0.0001,fld_Ke,np), ... % discretization for the lower branch of Phi_Vn
      linspace(fld_Ke,Kend,np)]'; % discretization for the upper branch of Phi_Vn
Phi_Ke = zeros(length(Ke),1);

% Compute function Phi_Ke on the 1st and 3rd branches of the curve Phi_Vn = 0
for i = 1:length(Ke)
    if i <= np
        Vn0cpBr = Vn0cpL; % take the lower branch of Vn0cp
    else
        Vn0cpBr = Vn0cpU; % take the upper branch of Vn0cp
    end
    Phi_Ke(i) = rhs_red_nan(Vn0cpBr(Ke(i)), Va0cp(Ke(i)), Ke(i), rhoN, ...
        rhoA, hp0, conc, 3); % @(x) -55
end

% Plot Phi_Ke restricted to lower and upper branches of the critical manifold
plot(Ke(1:np), Phi_Ke(1:np), '--', ...
    'Color', colors{1}, 'Linewidth', 2); % lower branch
plot(Ke(np:end), Phi_Ke(np:end), '--', ...
    'Color', colors{3}, 'Linewidth', 2); % piece of upper branch
plot([Ke(1) Ke(end)], [0 0], 'k--', 'Linewidth', 2);
%-------------------------------------------------------------------------%

% --------------------------- Zeros of Phi_Ke --------------------------- %
% Initial seeds to find the zeros of Phi_Ke
ind = find(Phi_Ke(1:end-1).*Phi_Ke(2:end) < 0); % indexes close to zeros
x0 = [Vn0cp(Ke(ind))'; Va0cp(Ke(ind))'; Ke(ind)']; % initial seeds close to zeros

% Compute zeros of Phi_Ke
for i = 1:size(x0,2)
    err = 1; % initial error
    iter = 1; % number of iterations
    while err >= 1e-12 && iter < 40
        [F, DF] = crtpts_sys(x0(:,i), rhoN, rhoA, hp0, conc); % System to find its eq. points
        x1 = x0(:,i) - DF\F; % Newton correction
        err = norm(x1 - x0(:,i))/norm(x1); % update error
        x0(:,i) = x1; % update seed
        iter = iter + 1;
    end
    err
    iter
    x0(:,i)
end

% Plot zeros of Phi_Ke restricted to the critical manifold
plot(x0(3,[1,3]), zeros(1, length(x0(3,[1,3]))), 'ko', 'MarkerFaceColor', 'k', ...
    'Linewidth', 2); % saddle points
plot(x0(3,2), 0, 'ko', 'MarkerFaceColor', 'w', 'Linewidth', 2); % unstable node/focus

% Create an inset figure
xstart = 0.48; xend = 0.855; ystart = 0.475; yend = 0.875;
axes('Position',[xstart ystart xend - xstart yend - ystart]);

hold on; box on; set(gca, 'Fontsize', 13);
xlim([5 25]); ylim([-0.01 0.06]);
plot(Ke1(1:np), Phi_Ke1(1:np), '-', 'Color', colors{1}, 'Linewidth', 2);
plot(Ke1(np:2*np+1), Phi_Ke1(np:2*np+1), '-', 'Color', colors{2}, 'Linewidth', 2);
plot(Ke1(2*np+1:end), Phi_Ke1(2*np+1:end), '-', 'Color', colors{3}, 'Linewidth', 2);
plot(Ke(1:np), Phi_Ke(1:np), '--', 'Color', colors{1}, 'Linewidth', 2);
plot(Ke(np:end), Phi_Ke(np:end), '--', 'Color', colors{3}, 'Linewidth', 2);
plot([Ke(1) Ke(end)], [0 0], 'k--', 'Linewidth', 2); % zero line
plot(x0(3,[1,3]), zeros(1, length(x0(3,[1,3]))), 'ko', 'MarkerFaceColor', 'k', ...
    'Linewidth', 2); % saddle points
plot(x0(3,2), 0, 'ko', 'MarkerFaceColor', 'w', 'Linewidth', 2); % unstable node/focus

% Save data
name_file = ['eqpoints_slow_subsystem_fldKe' '_rhon', num2str(rhoN), ...
             '_rhoa', num2str(rhoA), '.mat'];
save(name_file, 'x0', 'fld_Vn', 'fld_Ke', 'fld_Vn2', 'fld_Ke2');

% Save figure
name_fig = ['eqpoints_slow_subsystem_rhon', num2str(rhoN), ...
            '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Intersection points in the 3D implicit plot
plot3(axf3, Vn0cp(x0(3,:)), Va0cp(x0(3,:)), x0(3,:), 'ko', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'k');

% ------------------ Eigenvalues of the slow subsystem ------------------ %
pm = [1, -1];
ylabels = {'\lambda_+', '\lambda_-'};
lnstyles = {'-', '--'};

figure; hold on; box on; grid on; set(gca, 'Fontsize', 13); ax = gca;
set(get(gca, 'YAxis'), 'Exponent', -1); % set exponents exps = [0, -2];
xlabel('c', 'FontSize', 15, 'FontWeight', 'bold'); xlim([0 0.1]);
ylabel([ylabels{1}, ', ', ylabels{2}], 'Fontsize', 15, 'Fontweight', 'bold', 'Color', 'k');
plot([0 0.1], [0 0], ':', 'Color', [0 0 0], 'Linewidth', 1.5); % zero line
title('Eigenvalues of the slow subsystem', 'Fontsize', 14, ...
      'Fontweight', 'bold');

% col = {[0.8, 0.2 0.5], [0.2, 0.8 0.5], [0.5, 0.2 0.8]};
cmap = hsv(100); col = {cmap(30,:), cmap(50,:), cmap(10,:)};
pl = cell(3,2);
for l = 1:2
    for i = 1:3
        % Implicit differentiation of Vn([K^+]_e)
        [~, dVPhi_Vn, dKPhi_Vn] = diff_PhiVn(x0([1,3],i)', rhoN, hp0, conc);
        dm0dKe = -dKPhi_Vn./dVPhi_Vn;

        % Implicit differentiation of Va([K^+]_e)
        [~, dVPhi_Va, dKPhi_Va] = diff_PhiVa(x0([2,3],i)', rhoA, conc);
        dn0dKe = -dKPhi_Va./dVPhi_Va;

        % Derivatives of Phi_Ke
        [~, dVnPhi_Ke, dVaPhi_Ke, dKePhi_Ke] = diff_PhiKe(x0(:,i)', rhoN, rhoA, conc);

        % Compute derivative of function H(z) = h(x(z), y(z), z)
        dHdz = dVnPhi_Ke*dm0dKe + dVaPhi_Ke*dn0dKe + dKePhi_Ke;

        % Compute eigenvalues of the slow subsystem as a function of c
        cc = linspace(0, 0.1, 1500); D = 1;
        lam = 0.5*(cc/D + pm(l)*sqrt((cc/D).^2 - 4/D*dHdz));

        pl{i,l} = plot(cc, lam, 'Color', col{i}, 'Linestyle', lnstyles{l}, 'Linewidth', 2);
    end
end

% First legend
leg = legend(ax, [pl{:,1}], {'p_{l,1}', 'p_{l,2}', 'p_r'}, 'Fontsize', 14, ...
            'Location', 'NorthWest', 'NumColumns', 1);
leg.Title.String = ylabels{1};
leg.Title.FontSize = 14;
leg.ItemTokenSize = [20, 18];

% Superpose second figure
ax2 = axes('Position', get(ax, 'Position'), 'Visible', 'off');

% Second legend
leg2 = legend(ax2, [pl{:,2}], {'p_{l,1}', 'p_{l,2}', 'p_r'}, 'Fontsize', 14, ...
              'Location', 'NorthEast', 'NumColumns', 1);
leg2.Title.String = ylabels{2};
leg2.Title.FontSize = 14;
leg2.ItemTokenSize = [20, 18];
leg2.Position = leg.Position + [0.113 0 0 0];

% Save figure
name_fig = ['eigenvalues_slow_subsystem_rhoN', num2str(rhoN), ...
            '_rhoA', num2str(rhoA), '_D', num2str(D), '.eps'];
print(gcf, '-depsc', '-tiff', name_fig);
%-------------------------------------------------------------------------%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase portrait %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: wave speed and diffusion coefficient
c = 0.05; D = 1;

%------- Stability of the equilibrium points of the slow subsystem -------%
eqp = x0';

% Implicit differentiation of Vn([K^+]_e)
[~, dVnPhi_Vn, dKePhi_Vn] = diff_PhiVn(eqp(:,[1,3]), rhoN, hp0, conc);
dVndKe = -dKePhi_Vn./dVnPhi_Vn;

% Implicit differentiation of Va([K^+]_e)
[~, dVaPhi_Va, dKePhi_Va] = diff_PhiVa(eqp(:,[2,3]), rhoA, conc);
dVadKe = -dKePhi_Va./dVaPhi_Va;

[~, dVnPhi_Ke, dVaPhi_Ke, dKePhi_Ke] = diff_PhiKe(eqp, rhoN, rhoA, conc);
dPhi_Ke = dVnPhi_Ke.*dVndKe + dVaPhi_Ke.*dVadKe + dKePhi_Ke

[cc, dd] = meshgrid(-10:0.1:10, 1:0.1:30); cc1 = -10:0.1:10;
figure; hold on; box on; set(gca, 'Fontsize', 13);
lam1 = 0.5*(cc./dd + sqrt((cc./dd).^2 - 4*dPhi_Ke(3)./dd));
lam1(imag(lam1) ~= 0) = nan;
surf(cc, dd, lam1, 'EdgeColor', 'none');
lamp = 0.5*(cc1 + sqrt(cc1.^2 - 4*dHdz));
plot3(cc1, ones(length(cc1),1), lamp, 'k-', 'Linewidth', 1.5);
title('lam_+', 'Fontsize', 15, 'Fontweight', 'bold');

figure; hold on; box on; set(gca, 'Fontsize', 13);
lam2 = 0.5*(cc./dd - sqrt((cc./dd).^2 - 4*dPhi_Ke(3)./dd));
lam2(imag(lam2) ~= 0) = nan;
surf(cc, dd, lam2, 'EdgeColor', 'none');
lamn = 0.5*(cc1 - sqrt(cc1.^2 - 4*dHdz));
plot3(cc1, ones(length(cc1),1), lamn, 'k-', 'Linewidth', 1.5);
title('lam_-', 'Fontsize', 15, 'Fontweight', 'bold');
%-------------------------------------------------------------------------%

% Phase portrait restricted to [5 20] x [-0.1 0.1]
figure; hold on; box on; set(gca, 'Fontsize', 13);
xlabel('z', 'FontSize', 15, 'FontWeight', 'bold'); xlim([5 20]);
ylabel('w', 'FontSize', 15, 'FontWeight', 'bold'); ylim([-0.1 0.1]);

% Function Phi_Ke restricted to the critical manifold
h = @(z) rhs_red_nan(Vn0cp(z), Va0cp(z), z, rhoN, rhoA, hp0, conc, 3);

% ODE solver properties
Ps = @(t, x) eventPs(t, x, [1, fld_Ke, 215, 0.5, -0.5], [1, 1, 1, 2, 2]);
options_ode = odeset('Events', Ps); % 'AbsTol', 1e-5, 'RelTol', 1e-3, 

% Plotting streamlines of the slow subsystem: Region[5 20] x [-0.1 0.1]
[v, u] = meshgrid(5:2:fld_Ke, -0.1:0.01:0.1); % mesh of initial conditions
vp = u; up = (1/D)*(c*u - h(v)); % vector field of slow subsystem
% quiver(v, u, vp, up, 'Linewidth', 1.5, 'Color', [0.7 0.7 0.7]);
verts = stream2(v, u, vp, up, v, u);
% streamline(verts);
% [verts, ~] = streamslice(v, u, vp, up, 'arrows'); % vertices of streamlines
for i = 1:length(verts)
    ini = verts{i};
    [~, xst] = ode15s(@(t,x) [x(2); (1/D)*(c*x(2) - h(x(1)))], ...
        [0 2000], ini(1,:), options_ode);
    plot(xst(:,1), xst(:,2), 'Linewidth', 1);
    [~, xst] = ode15s(@(t,x) [x(2); (1/D)*(c*x(2) - h(x(1)))], ...
        [0 -2000], ini(1,:), options_ode);
    plot(xst(:,1), xst(:,2), 'Linewidth', 1);
    drawnow;
end

% Plotting streamlines of the slow subsystem: Region [195 215] x [-0.1 0.1]
[v, u] = meshgrid(195:2:215, -0.1:0.01:0.1); % mesh of initial conditions
vp = u; up = (1/D)*(c*u - h(v)); % vector field of slow subsystem
verts = stream2(v, u, vp, up, v, u);
for i = 1:length(verts)
    ini = verts{i};
    [tst, xst] = ode15s(@(t,x) [x(2); (1/D)*(c*x(2) - h(x(1)))], ...
        [0 2000], ini(1,:), options_ode);
    plot(xst(:,1), xst(:,2), 'Linewidth', 1);
    [~, xst] = ode15s(@(t,x) [x(2); (1/D)*(c*x(2) - h(x(1)))], ...
        [0 -2000], ini(1,:), options_ode);
    plot(xst(:,1), xst(:,2), 'Linewidth', 1);
    drawnow;
end

% Plotting equilibrium points
plot(x0(3,[1,3]), zeros(1, length(x0(3,[1,3]))), 'ko', 'Linewidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'k'); % saddle points
plot(x0(3,2), 0, 'ko', 'Linewidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'w'); % unstable node/focus

% Plotting folds
plot([fld_Ke fld_Ke], [-5, 5], 'k:', 'Linewidth', 2); % 1st asymptote (fold)
plot([fld_Ke2 fld_Ke2], [-5, 5], 'k:', 'Linewidth', 2); % 2nd asymptote (fold)

% Plotting invariant manifolds of both saddle points
[~, xu, xs] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, 1, 1, ...
                                fld_Ke, c, D, parm);
plot(xu(:,1), xu(:,2), 'k--', 'Linewidth', 2);
plot(xs(:,1), xs(:,2), 'k-', 'Linewidth', 2);
[~, xu, xs] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, -1, 1, ...
                                fld_Ke, c, D, parm);
plot(xu(:,1), xu(:,2), 'k--', 'Linewidth', 2);
plot(xs(:,1), xs(:,2), 'k-', 'Linewidth', 2);

[~, xu, xs] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, 1, -1, ...
                                fld_Ke, c, D, parm);
plot(xu(:,1), xu(:,2), 'k-', 'Linewidth', 2);
plot(xs(:,1), xs(:,2), 'k--', 'Linewidth', 2);
[~, xu, xs] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, -1, -1, ...
                                fld_Ke, c, D, parm);
plot(xu(:,1), xu(:,2), 'k-', 'Linewidth', 2);
plot(xs(:,1), xs(:,2), 'k--', 'Linewidth', 2);
%%
%-------- Illustration of the computation of a heteroclinic orbit --------%
% Figure of stable/unstable manifolds of saddle points
figure;
ax1 = subplot(1,2,1);
hold on; set(gca, 'FontSize', 13);
xlabel('z (lower branch)', 'FontSize', 15, 'FontWeight', 'bold'); xlim([10 fld_Ke]);
ylabel('w', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 1.25]);

ax2 = subplot(1,2,2);
hold on; set(gca, 'FontSize', 13, 'YTick', [], 'YColor', 'none');
xlabel('z (upper branch)', 'FontSize', 15, 'FontWeight', 'bold'); xlim([fld_Ke 215]);
linkaxes([ax1, ax2], 'y'); % synchronize limits

pos1 = get(ax1, 'Position'); set(ax1, 'Position', pos1 + [0 0 0.05254 0]);
pos2 = get(ax2, 'Position'); set(ax2, 'Position', pos2 + [-0.05254 0 0.05254 0]);

text(fld_Ke, 0.98, 'z = z^R', 'Fontsize', 16, 'Rotation', 90, ...
    'VerticalAlignment', 'bottom'); % fold annotation
xline(ax2, fld_Ke, 'k--', 'Linewidth', 1.5); % fold in the lower branch

% Plotting equilibrium points
% --> Lower branch
plot(ax1, x0(3,[1,3]), zeros(1, length(x0(3,[1,3]))), 'ko', 'MarkerFaceColor', 'k', ...
    'Linewidth', 1.5); % saddle points
plot(ax1, x0(3,2), 0, 'ko', 'MarkerFaceColor', 'w', 'Linewidth', 1.5); % unstable node/focus

% --> Upper branch
plot(ax2, x0(3,[1,3]), zeros(1, length(x0(3,[1,3]))), 'ko', 'MarkerFaceColor', 'k', ...
    'Linewidth', 1.5); % saddle points
plot(ax2, x0(3,2), 0, 'ko', 'MarkerFaceColor', 'w', 'Linewidth', 1.5); % unstable node/focus

cInt = zeros(2,1); % interval of values of c for which there is a heteroclinic orbit
frs_sgn = true; % to store the first value of c
chn_sgn = false; % to detect a change of sign

cc = [0.04 0.09]; % range of values c for which there is a heteroclinic
al = linspace(0.2,1,length(cc)); % alpha transparency
for i = 1:length(cc)
    c = cc(i)
    
    % Positive branches (c > 0) of the stable and unstable manifolds of 
    % the 1st and 2nd saddle points
    [~, Wu, Ws] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, ...
                                    sign(c), -sign(c), fld_Ke, c, D, parm);
    
    % If both Wu and Ws intersect with the Poincare section
    if ~chn_sgn
        w_Wu = Wu(end,2); % second component of the last integration point of Wu
        w_Ws = Ws(end,2); % second component of the last integration point of Ws
        if frs_sgn
            frs_sgn = false;
            sg = sign(w_Wu - w_Ws); % store the sign of the function to find a zero
        end
        chn_sgn = sign(w_Wu - w_Ws) ~= sg; % there is a change of sign (root assured)
        if chn_sgn
            cInt(2) = c; % store the second value of c
            
            % Plotting invariant manifolds (after change of sign)
            plot(ax1, Ws(:,1), Ws(:,2), '-', 'Color', [0 0.447 0.741 1], ...
                'Linewidth', 1.5);
            plot(ax2, Wu(:,1), Wu(:,2), '-', 'Color', [0.850 0.325 0.098 1], ...
                'Linewidth', 1.5);
            
            % Plotting intersection points (after change of sign)
            scatter(ax1, Ws(end,1), Ws(end,2), 'MarkerFaceColor', [0 0.447 0.741], ...
                'MarkerEdgeColor', [0 0.447 0.741], 'MarkerFaceAlpha', 1, ...
                'Linewidth', 1.5);
            scatter(ax2, Wu(end,1), Wu(end,2), 'MarkerFaceColor', [0.850 0.325 0.098], ...
                'MarkerEdgeColor', [0.850 0.325 0.098], 'MarkerFaceAlpha', 1, ...
                'Linewidth', 1.5);
        else
            cInt(1) = c; % update the first value of c
            
            % Plotting invariant manifolds (before change of sign)
            if exist('pl1', 'var')
                delete(pl1); delete(pl2);
                delete(ps1); delete(ps2);
            end
            pl1 = plot(ax1, Ws(:,1), Ws(:,2), '-', 'Color', [0 0.447 0.741 0.3], ...
                'Linewidth', 1.5);
            pl2 = plot(ax2, Wu(:,1), Wu(:,2), '-', 'Color', [0.850 0.325 0.098 0.3], ...
                'Linewidth', 1.5);
            
            % Plotting intersection points (before change of sign)
            ps1 = scatter(ax1, Ws(end,1), Ws(end,2), 'MarkerFaceColor', [0 0.447 0.741], ...
                'MarkerEdgeColor', [0 0.447 0.741], 'MarkerFaceAlpha', 0.3, ...
                'Linewidth', 1.5);
            ps2 = scatter(ax2, Wu(end,1), Wu(end,2), 'MarkerFaceColor', [0.850 0.325 0.098], ...
                'MarkerEdgeColor', [0.850 0.325 0.098], 'MarkerFaceAlpha', 0.3, ...
                'Linewidth', 1.5);
        end
    end
    
    drawnow;
end

% Finding a heteroclinic orbit by solving w_Wu - w_Ws = 0 with given tol
F = @(c) inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, ...
                           1, -1, fld_Ke, c, D, parm);
                       
% --> Bisection method
% c_ht = bisection_method(F, cInt(1), cInt(2), 5e-12)

% --> Secant method
c_ht = secant_method(F, cInt(1), cInt(2), 5e-12)

% --> Newton method
% figure; hold on; box on;
% argin = {x0(:,1)', x0(:,3)', Vn0cp, Va0cp, fld_Ke, D, parm};
% c_ht = newton_method(cInt(1), 5e-12, argin)

% Plotting singular heteroclinic orbit
[~, Wu, Ws] = inv_manifolds_slw(x0(:,1)', x0(:,3)', Vn0cp, Va0cp, ...
                        sign(c_ht), -sign(c_ht), fld_Ke, c_ht, D, parm);
plot(ax1, Ws(:,1), Ws(:,2), 'k-', 'Linewidth', 1.5);
plot(ax2, Wu(:,1), Wu(:,2), 'k-', 'Linewidth', 1.5);
plot(ax2, Wu(end,1), Wu(end,2), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 1.5);
sgt = sgtitle('Singular heteroclinic orbit', 'Fontsize', 15, 'Fontweight', 'bold');

% Save figure
name_fig = ['singular_heteroclinic_D', num2str(D), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '.eps'];
print(gcf, '-depsc', name_fig);

% Save singular heteroclinic orbit
name_file = ['singular_heteroclinic_D', num2str(D), '.txt'];
fclose(fopen(name_file, 'w+'));

file = fopen(name_file, 'a');
for i = 1:length(Ws) % save unstable manifold in the lower branch
    fprintf(file, '%16.15f %16.15f %16.15f %16.15f\r\n', [Vn0cpL(Ws(i,1)) ...
        Va0cp(Ws(i,1)) Ws(i,:)]);
end

for i = length(Wu):-1:1 % save stable manifold in the upper branch
    fprintf(file, '%16.15f %16.15f %16.15f %16.15f\r\n', [Vn0cpU(Wu(i,1)) ...
        Va0cp(Wu(i,1)) Wu(i,:)]);
end
fclose(file);
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

eps = 0;

% Initial conditions for the neuron-astrocyte network model
x0 = 3.5*ones(N,1);

% Integrating neuron-astrocyte network model
options_ode = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
pumpK = true; % reset pumpK variable again to true
t0 = 0; tf = 25000; % initial/final time of the integration
[t, y] = ode15s(@(t, x) red_bis_nan_model(t, x, N, rhoN, rhoA, hp0, conc, ...
                            {Vn0cpL, Vn0cpU}, Va0cp, fld_Ke), [t0 tf], x0, options_ode);

% [K^+]_e(x,t) vs t
figure; hold on; grid on; box on; set(gca, 'FontSize', 13);
xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); 
ylabel('[K^+]_e (mM)', 'FontSize', 15, 'FontWeight', 'bold');
title('[K^+]_e time evolution', 'FontSize', 15, 'FontWeight', 'bold');
h = plot(t/1000, y(:,1:N), 'Linewidth', 1.5);
set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines

% -- Colorbar
c = colorbar; colormap(parula(N)); c.FontSize = 13;
c.TickLabels = cellfun(@num2str, num2cell(0:10:N), 'UniformOutput', false); 
set(c.Label, 'String', '# cell', 'Fontsize', 15, 'Fontweight', 'bold');

% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Ke_time.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Neurons' extracellular K^+ concentration
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Extracellular K^+ concentration', 'Fontsize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, y(:,1:N)); shading interp; caxis([5 150]);
% figure; imagesc([1, N], [t(1), t(end)], flipud(y(:,1:N))); shading interp; caxis([3.5 70]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb.Ruler, 'Exponent', 1); % set exponent 1 to the colorbar
set(cb.Ruler.SecondaryLabel, 'Position', [1 154 0]); % location of exponent
set(cb, 'Ticks', [5 150]); % skip intermediate ticks on colorbar
set(cb.Label, 'String', '[K^+]_e (mM)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar

% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Ke.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

%%%%%%%%%%%%%%%%%% Plotting Vn and Va in terms of [K^+]_e %%%%%%%%%%%%%%%%%
V_N = zeros(length(t), N); V_A = zeros(length(t), N);

for i = 1:N
    V_N(:,i) = Vn0cp(y(:,i));
    V_A(:,i) = Va0cp(y(:,i));
end

% Neurons' membrane potential (as function of [K^+]_e)
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Neurons'' membrane potential', 'FontSize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, V_N(:,1:N)); shading interp; caxis([-70 20]);
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', 'V_N (mV)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar

% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Vn.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Vn(x,t) vs t
figure; hold on; grid on; box on; set(gca, 'FontSize', 13);
xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); 
ylabel('V_N (mV)', 'FontSize', 15, 'FontWeight', 'bold');
title('V_N time evolution', 'FontSize', 15, 'FontWeight', 'bold');
h = plot(t/1000, V_N(:,1:N), 'Linewidth', 1.5);
set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines

% -- Colorbar
c = colorbar; colormap(parula(N)); c.FontSize = 13;
c.TickLabels = cellfun(@num2str, num2cell(0:10:N), 'UniformOutput', false); 
set(c.Label, 'String', '# cell', 'Fontsize', 15, 'Fontweight', 'bold');

% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Vn_time.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Astrocytes' membrane potential (as function of [K^+]_e)
figure; hold on; box on; set(gca, 'FontSize', 13);
xlabel('# cell', 'FontSize', 15, 'FontWeight', 'bold'); xlim([1 N]);
ylabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); ylim([0 tf/1000]);
title('Astrocytes'' membrane potential', 'FontSize', 15, 'FontWeight', 'bold');
surf(1:N, t/1000, V_A(:,1:N)); shading interp; caxis([-80 0]); % -90 20
set(gca, 'Layer', 'top'); % placement of grid lines/tick marks
axpos = get(gca, 'Position'); % get current position
set(gca, 'Position', axpos + [-0.035 0 0 0]);
cb = colorbar; % add colorbar to the Poincare section plot
set(cb, 'Position', cb.Position + [0.1 0 0 0]);
set(cb, 'Ticks', cb.Ticks([1,end])); % skip intermediate ticks on colorbar
set(cb.Label, 'String', 'V_A (mV)', 'Fontsize', 14, 'Fontweight', 'bold', ...
    'VerticalAlignment', 'baseline'); % set label to the colorbar
    
% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Va.eps'];
print(gcf, '-depsc', '-tiff', name_fig);

% Va(x,t) vs t
figure; hold on; grid on; box on; set(gca, 'FontSize', 13);
xlabel('t (s)', 'FontSize', 15, 'FontWeight', 'bold'); 
ylabel('V_A (mV)', 'FontSize', 15, 'FontWeight', 'bold');
title('V_A time evolution', 'FontSize', 15, 'FontWeight', 'bold');
h = plot(t/1000, V_A(:,1:N), 'Linewidth', 1.5);
set(h, {'Color'}, num2cell(parula(N),2)); % change default color of lines

% -- Colorbar
c = colorbar; colormap(parula(N)); c.FontSize = 13;
c.TickLabels = cellfun(@num2str, num2cell(0:10:N), 'UniformOutput', false); 
set(c.Label, 'String', '# cell', 'Fontsize', 15, 'Fontweight', 'bold');

% Save figure
name_fig = ['csd_red_nan_eps', num2str(eps), '_rhon', num2str(rhoN), ...
    '_rhoa', num2str(rhoA), '_Va_time.eps'];
print(gcf, '-depsc', '-tiff', name_fig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%