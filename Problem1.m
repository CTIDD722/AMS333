%% Problem 1 Yeast growth
% author: Yunpeng Chu
% Date:   2025/9/17
% Note:   Please use version 2013a or above

clear; clc; close all;
%% Common parameters
tau  = 1.5;                 % doubling time 
rate = log(2)/tau;          % growth rate 

%% (a)
N0_a    = 1;                        % initial cell number
t_final = 72;                       % total time
dt_a    = 3;                        % sampling interval
tgrid_A  = 0:dt_a:t_final;          % sampling times
N_expA  = N0_a * 2.^(tgrid_A./tau); % unrestricted growth model
N_3days = N_expA(end);

fprintf('===== (a) =====\n\n');
fprintf('Estimated number of cells after 3 days N = %.4e\n\n', N_3days);

figure('Name','(a) Exponential growth','Color','w');
semilogy(tgrid_A, N_expA, 'o-','LineWidth',1.5); grid on;
xlabel('Time t (hours)'); ylabel('Cell number N(t)'); title('Unrestricted cell growth (3 days)');

%% (b)
d_um     = 6;               % diameter 
r_m      = (d_um/2) * 1e-6; % radius 
Vcell_m3 = (4/3)*pi*r_m^3;  % m^3
Vcell_L  = Vcell_m3 / 1e-3; % L

V_total = N_3days * Vcell_L;

fprintf('===== (b) =====\n\n');
fprintf('Single cell volume  V_cell = %.3e L\n', Vcell_L);
fprintf('Total volume after 3 days  V_total = %.3e L\n\n', V_total);

%% (c)
N0_c  = 1e4; % initial cell number
V_max = 1.0; % maximum volume
N_max = V_max / Vcell_L;    

t_max = tau * log2(N_max / N0_c);
fprintf('===== (c) =====\n\n');
fprintf('Estimated by "space limit": N_max ≈ %.3e cells\n', N_max);
fprintf('Time for exponential model to reach this limit t ≈ %.2f hours (≈ %.2f days)\n\n', t_max, t_max/24);

%% (d) 
dens_max = 2e8;                     % cells/mL
K        = dens_max * 1000 * V_max; % carrying capacity for 1 L
V_K_max = K * Vcell_L;       
fraction = V_K_max / V_max;

fprintf('===== (d) =====\n\n');
fprintf('Carrying capacity for 1 L K = %.2e cells\n', K);
fprintf('Total cell volume at maximum density ≈ %.3f mL (≈ %.2f%% of 1 L)\n\n', ...
        V_K_max*1e3, fraction*100);

%% (e)
N0_e   = N0_c;
tmax_e = max(t_max*1.4, 100);         
tgrid_E = linspace(0, tmax_e, 600);

A        = (K - N0_e)/N0_e;
N_logE   = K ./ (1 + A .* exp(-rate.*tgrid_E));

% Plot
figure('Name','(e) Logistic model','Color','w');
semilogy(tgrid_E, N_logE, 'LineWidth',1.8); grid on;
xlabel('Time t (hours)'); ylabel('Cell number N(t)');
title('(e) Logistic growth: N_0=10^4, K=2\times10^{11}');
yline(K,'--','K');

% Find 0.99K point
alpha    = 0.99;
t_99     = (1/rate) * log( (alpha/(1-alpha)) * (K/N0_e - 1) );
days_99K = t_99/24;

figure('Name','(e) Logistic model','Color','w'); hold on; grid on; box on;
semilogy(tgrid_E, N_logE, 'LineWidth',1.8, 'Color',[0 0.45 0.74]);
yline(K,'--k','K');
patch([0 tmax_e tmax_e 0], [alpha*K alpha*K K K], ...
      [1 0 0], 'FaceAlpha',0.08, 'EdgeColor','none');

xline(t_99,'-r', ...
      'Label',sprintf('0.99K at %.1f h',t_99), ...
      'LabelOrientation','horizontal', 'LineWidth',1.6);
semilogy(t_99, alpha*K, 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
xlabel('Time t (hours)'); 
ylabel('Cell number N(t)');
title('(e) Logistic growth: N_0=10^4, K=2\times10^{11}');
xlim([0 tmax_e]);
ylim([min(N_logE)/1.5, K*1.2]);
set(gca,'YMinorGrid','on');

fprintf('===== (e) =====\n\n');
fprintf('Logistic: t_{0.99K} ≈ %.2f hours (≈ %.2f days)\n', t_99, days_99K);


