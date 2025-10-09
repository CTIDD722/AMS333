%% Project 2
%  Author:Yunpeng Chu
%  Date:  2025/10/09
%  Note:  please make sure the version of MATLAB is above R2016a

clear; clc; close all;
%% Problem 2
n_birth = 18;   % number of litters per year
S_1st   = 0.33; % survival rate in the first month
S_conti = 0.95; % monthly survival rate after the first month
Rf      = 0.5;  % female ratio
Month   = 11;   % months 2–12

S_adult   = S_conti^12;
S_2to12   = S_conti^Month;
R_recruit = Rf * n_birth * S_1st * S_2to12;
D_adult   = 1-S_adult;
Ra        = S_adult + R_recruit;  % annual net growth factor
alpha     = log(Ra);              % annual per-capita growth rate
fprintf('\n[Problem 2] alpha = %.9f\n', alpha);
fprintf('[Problem 2] D_adult = %.9f\n', D_adult);

%% Problem 3
S1   = 0.7;        % monthly survival without prey
S_yr = S1^12;
D_yr = 1-S_yr;
beta = -log(S_yr); % annual per-capita death rate
fprintf('[Problem 3] beta = %.9f\n', beta);
fprintf('[Problem 3] D_yr = %.9f\n', D_yr);

%% Problem 4
gamma = 365 / 1000;        % per-year per-density predation rate
epsi  = (0.10 * 1.5) / 10; % conversion efficiency
fprintf('[Problem 4] gamma = %.9f, epsilon = %.9f\n', gamma, epsi);

%% Problem 5
dt   = 0.001; % time step
Tend = 40;    % total time
t    = 0:dt:Tend;

% equilibrium
P_eq = beta / (epsi * gamma);
Z_eq = alpha / gamma;
fprintf('[INFO] Equilibrium (P*,Z*) = (%.9f, %.9f)\n', P_eq, Z_eq);

%% Problem 6-7
% initial conditions
cases = { ...
    struct('name',"Case A (baseline)", 'P0',400, 'Z0',1), ...
    struct('name',"Case B (800, 2)",  'P0',800, 'Z0',2), ...
    struct('name',"Case C (200, 0.5)",'P0',200, 'Z0',0.5) ...
};
results = cell(1, numel(cases));

for ci = 1:numel(cases)
    C = cases{ci};
    [P, Z] = simulate(C.P0, C.Z0, alpha, beta, gamma, epsi, t, dt);
    M = compute_M(P, Z, t, alpha, beta, gamma, epsi);
    results{ci} = struct('name',C.name,'P0',C.P0,'Z0',C.Z0,'P',P,'Z',Z,'M',M);
    timeseries(t, P, Z, C.name);                           
    plot_phase_plane(P, Z, alpha, beta, gamma, epsi, P_eq, Z_eq, C.name); 
end

fprintf('\n===== Problem 7: Similarities & Differences =====\n');
R_A = results{1}; R_B = results{2}; R_C = results{3};

% similarities
fprintf('Similarities:\n');
fprintf('  • All three cases produce closed orbits around the same equilibrium (P*,Z*),\n');
fprintf('    with predator peaks lagging prey peaks by ~%.2f–%.2f years.\n', ...
    min([R_A.M.lag_mean,R_B.M.lag_mean,R_C.M.lag_mean]), ...
    max([R_A.M.lag_mean,R_B.M.lag_mean,R_C.M.lag_mean]));

% differences
fprintf('Differences (relative to baseline Case A):\n');
print_delta('Prey period (yr)', R_A.M.PerP, R_B.M.PerP, R_C.M.PerP);
print_delta('Pred period (yr)', R_A.M.PerZ, R_B.M.PerZ, R_C.M.PerZ);
print_delta('Prey mean amplitude', R_A.M.AmpP, R_B.M.AmpP, R_C.M.AmpP);
print_delta('Pred mean amplitude', R_A.M.AmpZ, R_B.M.AmpZ, R_C.M.AmpZ);
print_delta('Time P<1/100 (yr)', R_A.M.time_below_P, R_B.M.time_below_P, R_C.M.time_below_P);
print_delta('Time Z<1/100 (yr)', R_A.M.time_below_Z, R_B.M.time_below_Z, R_C.M.time_below_Z);
fprintf('  • Max growth/decline rates scale with amplitude; see per-case printouts above.\n');

%% Problem 8: Logistic prey growth
K = 3000;  % carrying capacity of hares
fprintf('\n===== Problem 8: Logistic =====\n', K);

% Logistic 
P_eq_log = P_eq;
Z_eq_log = (alpha/gamma) * (1 - P_eq_log / K);
fprintf('[Equilibria] LV  (P*,Z*)=(%.4f, %.4f) | Logistic (P*,Z*)=(%.4f, %.4f)\n', ...
    P_eq, Z_eq, P_eq_log, Z_eq_log);

for ci = 1:numel(cases)
    C = cases{ci};
    P_lv = results{ci}.P;  Z_lv = results{ci}.Z;  M_lv = results{ci}.M;
    [P_log, Z_log] = simulate_log(C.P0, C.Z0, alpha, beta, gamma, epsi, K, t, dt);
    M_log = compute_M_log(P_log, Z_log, t, alpha, beta, gamma, epsi, K);
    % comparison plots
    timeseries_log(t, P_lv, Z_lv, P_log, Z_log, sprintf('%s — LV vs Logistic', C.name));
    plot_phase_plane_logistic(P_log, Z_log, alpha, beta, gamma, epsi, K, P_eq_log, Z_eq_log, sprintf('%s — Logistic', C.name));
    fprintf('\n[Problem 8 | %s] Logistic vs LV:\n', C.name);
    print_cmp('Prey period (yr)',    M_lv.PerP, M_log.PerP);
    print_cmp('Pred period (yr)',    M_lv.PerZ, M_log.PerZ);
    print_cmp('Prey mean amplitude', M_lv.AmpP, M_log.AmpP);
    print_cmp('Pred mean amplitude', M_lv.AmpZ, M_log.AmpZ);
    print_cmp('Prey min/max',        [M_lv.P_min M_lv.P_max], [M_log.P_min M_log.P_max]);
    print_cmp('Pred min/max',        [M_lv.Z_min M_lv.Z_max], [M_log.Z_min M_log.Z_max]);
    print_cmp('Time P<1/100 (yr)',   M_lv.time_below_P, M_log.time_below_P);
    print_cmp('Time Z<1/100 (yr)',   M_lv.time_below_Z, M_log.time_below_Z);
    if ~isnan(M_lv.lag_mean) && ~isnan(M_log.lag_mean)
        fprintf('  Predator peak lag (mean): LV=%.3f yr | Logistic=%.3f yr\n', M_lv.lag_mean, M_log.lag_mean);
    else
        fprintf('  Predator peak lag: insufficient peaks in one model.\n');
    end
    if (isnan(M_log.PerP) || M_log.AmpP < 0.5*M_lv.AmpP) || (isnan(M_log.PerZ) || M_log.AmpZ < 0.5*M_lv.AmpZ)
        fprintf('  Obs: Logistic prey reduces oscillation amplitude and often damps to coexistence equilibrium.\n');
    else
        fprintf('  Obs: Logistic prey retains cycles but with altered amplitude/period vs LV.\n');
    end
end

%% Problem 9
kappa_h = 4; % handling time
kappa_y = kappa_h / (24*365);     
fprintf('\n===== Problem 9: Holling II =====\n', kappa_y);

% existence check
den_h = (epsi - beta*kappa_y);
if den_h <= 0
    fprintf('Warning: epsi - beta*kappa <= 0, no positive coexistence equilibrium; predators cannot persist.\n');
end

for ci = 1:numel(cases)
    C   = cases{ci};
    Plv = results{ci}.P;  Zlv = results{ci}.Z;  Mlv = results{ci}.M;
    [Phol, Zhol] = simulate_holl(C.P0, C.Z0, alpha, beta, gamma, epsi, kappa_y, t, dt);
    Mhol = compute_M_holl(Phol, Zhol, t, alpha, beta, gamma, epsi, kappa_y);
    % comparison plots
    timeseries_holling(t, Plv, Zlv, Phol, Zhol, sprintf('%s — LV vs Holling-II', C.name));
    if den_h > 0
        P_eq = beta / (gamma*den_h);
        Z_eq = (alpha/gamma) * (1 + gamma*kappa_y*P_eq);
    else
        P_eq = NaN; Z_eq = NaN;
    end
    phaseplane_holling(Phol, Zhol, alpha, beta, gamma, epsi, kappa_y, P_eq, Z_eq, sprintf('%s — Holling-II', C.name));
    fprintf('\n[Problem 9 | %s] Holling-II vs LV:\n', C.name);
    print_name('Prey period (yr)',    Mlv.PerP, Mhol.PerP, 'Hol');
    print_name('Pred period (yr)',    Mlv.PerZ, Mhol.PerZ, 'Hol');
    print_name('Prey mean amplitude', Mlv.AmpP, Mhol.AmpP, 'Hol');
    print_name('Pred mean amplitude', Mlv.AmpZ, Mhol.AmpZ, 'Hol');
    print_name('Prey min/max',        [Mlv.P_min Mlv.P_max], [Mhol.P_min Mhol.P_max], 'Hol');
    print_name('Pred min/max',        [Mlv.Z_min Mlv.Z_max], [Mhol.Z_min Mhol.Z_max], 'Hol');
    print_name('Time P<1/100 (yr)',   Mlv.time_below_P, Mhol.time_below_P, 'Hol');
    print_name('Time Z<1/100 (yr)',   Mlv.time_below_Z, Mhol.time_below_Z, 'Hol');
    if ~isnan(Mlv.lag_mean) && ~isnan(Mhol.lag_mean)
        fprintf('  Predator peak lag (mean): LV=%.3f yr | Holling=%.3f yr\n', Mlv.lag_mean, Mhol.lag_mean);
    else
        fprintf('  Predator peak lag: insufficient peaks in one model.\n');
    end
    if (isnan(Mhol.PerP) || Mhol.AmpP < 0.6*Mlv.AmpP) || (isnan(Mhol.PerZ) || Mhol.AmpZ < 0.6*Mlv.AmpZ)
        fprintf('  Obs: Type-II saturation reduces effective intake, usually suppressing amplitude and converging faster.\n');
    else
        fprintf('  Obs: Cycles remain, but period/amplitude differ systematically from LV.\n');
    end
end

%% Problem 10
fprintf('\n===== Problem 10: Logistic + Holling =====\n');       
den_lh   = (epsi - beta*kappa_y);
if den_lh <= 0
    fprintf('Warning: epsi - beta*kappa <= 0\n');
end

for ci = 1:numel(cases)
    C   = cases{ci};
    Plv = results{ci}.P;  Zlv = results{ci}.Z;  Mlv = results{ci}.M;
    [Plh, Zlh] = simulate_logistic_holling(C.P0, C.Z0, alpha, beta, gamma, epsi, K, kappa_y, t, dt);
    Mlh = compute_M_logH(Plh, Zlh, t, alpha, beta, gamma, epsi, K, kappa_y);
    timeseries_logH(t, Plv, Zlv, Plh, Zlh, sprintf('%s — LV vs Logistic+Holling-II', C.name));
    if den_lh > 0
        Ueq_lh = beta / (gamma*den_lh);                                      
        Veq_lh = (alpha/gamma) * (1 + gamma*kappa_y*Ueq_lh) * (1 - Ueq_lh/K);
    else
        Ueq_lh = NaN; Veq_lh = NaN;
    end
    plot_phase_plane_logH(Plh, Zlh, alpha, beta, gamma, epsi, K, kappa_y, Ueq_lh, Veq_lh, sprintf('%s — Logistic+Holling-II', C.name));

    fprintf('\n[Problem 10 | %s] (Log+Hol) vs LV:\n', C.name);
    print_name('Prey period (yr)',    Mlv.PerP, Mlh.PerP, 'Log+Hol');
    print_name('Pred period (yr)',    Mlv.PerZ, Mlh.PerZ, 'Log+Hol');
    print_name('Prey mean amplitude', Mlv.AmpP, Mlh.AmpP, 'Log+Hol');
    print_name('Pred mean amplitude', Mlv.AmpZ, Mlh.AmpZ, 'Log+Hol');
    print_name('Prey min/max',        [Mlv.P_min Mlv.P_max], [Mlh.P_min Mlh.P_max], 'Log+Hol');
    print_name('Pred min/max',        [Mlv.Z_min Mlv.Z_max], [Mlh.Z_min Mlh.Z_max], 'Log+Hol');
    print_name('Time P<1/100 (yr)',   Mlv.time_below_P, Mlh.time_below_P, 'Log+Hol');
    print_name('Time Z<1/100 (yr)',   Mlv.time_below_Z, Mlh.time_below_Z, 'Log+Hol');
    if ~isnan(Mlv.lag_mean) && ~isnan(Mlh.lag_mean)
        fprintf('  Predator peak lag (mean): LV=%.3f yr | Log+Hol=%.3f yr\n', Mlv.lag_mean, Mlh.lag_mean);
    else
        fprintf('  Predator peak lag: insufficient peaks in one model.\n');
    end
    if (isnan(Mlh.PerP) || Mlh.AmpP < 0.6*Mlv.AmpP) || (isnan(Mlh.PerZ) || Mlh.AmpZ < 0.6*Mlv.AmpZ)
        fprintf('  Obs: The dual limits (logistic + handling) cause noticeable damping and faster convergence to coexistence.\n');
    else
        fprintf('  Obs: Period and amplitude differ from LV.\n');
    end
end

%% functions
function [P, Z] = simulate(P0, Z0, alpha, beta, gamma, epsi, t, dt)
nT = numel(t);
P  = zeros(1,nT); 
Z  = zeros(1,nT);
P(1) = P0; Z(1) = Z0;
for k = 1:nT-1
    dP = alpha * P(k) - gamma * P(k) * Z(k);
    dZ = epsi * gamma * P(k) * Z(k) - beta * Z(k);
    P(k+1) = max(P(k) + dt * dP, 0);
    Z(k+1) = max(Z(k) + dt * dZ, 0);
end
end

function [P, Z] = simulate_log(P0, Z0, alpha, beta, gamma, epsi, K, t, dt)
nT = numel(t);
P  = zeros(1,nT); 
Z  = zeros(1,nT);
P(1) = P0; Z(1) = Z0;
for k = 1:nT-1
    dP = alpha * P(k) * (1 - P(k)/K) - gamma * P(k) * Z(k);
    dZ = epsi * gamma * P(k) * Z(k) - beta * Z(k);
    P(k+1) = max(P(k) + dt * dP, 0);
    Z(k+1) = max(Z(k) + dt * dZ, 0);
end
end

function [P, Z] = simulate_holl(P0, Z0, alpha, beta, gamma, epsi, kappa, t, dt)
nT = numel(t); P = zeros(1,nT); Z = zeros(1,nT); P(1)=P0; Z(1)=Z0;
for k=1:nT-1
    func = (gamma*P(k)) / (1 + gamma*kappa*P(k));  
    dP = alpha*P(k) - func*Z(k);
    dZ = epsi*func*Z(k) - beta*Z(k);
    P(k+1) = max(P(k) + dt*dP, 0);
    Z(k+1) = max(Z(k) + dt*dZ, 0);
end
end

function M = compute_M(P, Z, t, alpha, beta, gamma, epsi)
% continuous derivatives
dPdt = alpha.*P - gamma.*P.*Z;
dZdt = epsi.*gamma.*P.*Z - beta.*Z;

% global min/max
M.P_min = min(P); M.P_max = max(P);
M.Z_min = min(Z); M.Z_max = max(Z);

% peaks and period/amplitude
burn_in = 5; maskBI = t >= burn_in;
[tPkP, ~, iPkP] = detect_peaks(P(maskBI), t(maskBI), 1.0);
[tPkZ, ~, iPkZ] = detect_peaks(Z(maskBI), t(maskBI), 1.0);
M.PerP = mean(diff(tPkP));  
if isempty(M.PerP), M.PerP = NaN; end
M.PerZ = mean(diff(tPkZ));  
if isempty(M.PerZ), M.PerZ = NaN; end
M.AmpP = mean_peak(P(maskBI), iPkP);
M.AmpZ = mean_peak(Z(maskBI), iPkZ);

% max growth and decline
[M.max_dPdt, iPpos] = max(dPdt);   M.tPpos = t(iPpos);
[M.min_dPdt, iPneg] = min(dPdt);   M.tPneg = t(iPneg);
[M.max_dZdt, iZpos] = max(dZdt);   M.tZpos = t(iZpos);
[M.min_dZdt, iZneg] = min(dZdt);   M.tZneg = t(iZneg);

% peak lag
lag_list = [];
for i = 1:numel(tPkP)
    tz = tPkZ(find(tPkZ >= tPkP(i), 1, 'first'));
    if ~isempty(tz), lag_list(end+1) = tz - tPkP(i); 
    end 
end
M.lag_mean   = mean(lag_list);
M.lag_median = median(lag_list);

% time below 1/100 km^2
thr = 0.01;
dtloc = t(2)-t(1);
M.time_below_P = dtloc * sum(P < thr);
M.time_below_Z = dtloc * sum(Z < thr);

fprintf('\n[Metrics] %s\n', evalin('caller','C.name'));
fprintf('  Prey  min/max: [%g, %g], period ~ %.3f yr, mean amp ~ %.3f\n', M.P_min, M.P_max, M.PerP, M.AmpP);
fprintf('  Pred  min/max: [%g, %g], period ~ %.3f yr, mean amp ~ %.3f\n', M.Z_min, M.Z_max, M.PerZ, M.AmpZ);
fprintf('  Max dP/dt = %+g at t=%.3f;  Min dP/dt = %+g at t=%.3f\n', M.max_dPdt, M.tPpos, M.min_dPdt, M.tPneg);
fprintf('  Max dZ/dt = %+g at t=%.3f;  Min dZ/dt = %+g at t=%.3f\n', M.max_dZdt, M.tZpos, M.min_dZdt, M.tZneg);
if ~isempty(lag_list)
    fprintf('  Predator peak lag: mean %.3f yr (median %.3f yr)\n', M.lag_mean, M.lag_median);
else
    fprintf('  Predator peak lag: insufficient peaks\n');
end
fprintf('  Time below 1/100 km^2: Prey %.3f yr, Pred %.3f yr\n', M.time_below_P, M.time_below_Z);
end

function M = compute_M_log(P, Z, t, alpha, beta, gamma, epsi, K)
% continuous derivatives
dPdt = alpha.*P.*(1 - P./K) - gamma.*P.*Z;
dZdt = epsi.*gamma.*P.*Z - beta.*Z;

% global min/max
M.P_min = min(P); M.P_max = max(P);
M.Z_min = min(Z); M.Z_max = max(Z);

% max growth and decline
burn_in = 5; maskBI = t >= burn_in;
[tPkP, ~, iPkP] = detect_peaks(P(maskBI), t(maskBI), 1.0);
[tPkZ, ~, iPkZ] = detect_peaks(Z(maskBI), t(maskBI), 1.0);
M.PerP = mean(diff(tPkP));  
if isempty(M.PerP), M.PerP = NaN; end
M.PerZ = mean(diff(tPkZ));  
if isempty(M.PerZ), M.PerZ = NaN; end
M.AmpP = mean_peak(P(maskBI), iPkP);
M.AmpZ = mean_peak(Z(maskBI), iPkZ);

% peaks and period/amplitude
[M.max_dPdt, iPpos] = max(dPdt);   M.tPpos = t(iPpos);
[M.min_dPdt, iPneg] = min(dPdt);   M.tPneg = t(iPneg);
[M.max_dZdt, iZpos] = max(dZdt);   M.tZpos = t(iZpos);
[M.min_dZdt, iZneg] = min(dZdt);   M.tZneg = t(iZneg);

% peak lag
lag_list = [];
for i = 1:numel(tPkP)
    tz = tPkZ(find(tPkZ >= tPkP(i), 1, 'first'));
    if ~isempty(tz), lag_list(end+1) = tz - tPkP(i); end %#ok<AGROW>
end
M.lag_mean   = mean(lag_list);
M.lag_median = median(lag_list);

% time below 1/100 km^2
thr = 0.01; dtloc = t(2)-t(1);
M.time_below_P = dtloc * sum(P < thr);
M.time_below_Z = dtloc * sum(Z < thr);
end

function M = compute_M_holl(P, Z, t, alpha, beta, gamma, epsi, kappa)
func = (gamma.*P) ./ (1 + gamma*kappa.*P); 
% continuous derivatives
dPdt = alpha.*P - func.*Z;
dZdt = epsi.*func.*Z - beta.*Z;

% global min/max
M.P_min = min(P); M.P_max = max(P);
M.Z_min = min(Z); M.Z_max = max(Z);

% max growth and decline
burn_in = 5; mask = t>=burn_in;
[tPkP,~,iPkP] = detect_peaks(P(mask), t(mask), 1.0);
[tPkZ,~,iPkZ] = detect_peaks(Z(mask), t(mask), 1.0);
M.PerP = mean(diff(tPkP)); if isempty(M.PerP), M.PerP=NaN; end
M.PerZ = mean(diff(tPkZ)); if isempty(M.PerZ), M.PerZ=NaN; end
M.AmpP = mean_peak(P(mask), iPkP);
M.AmpZ = mean_peak(Z(mask), iPkZ);

% peaks and period/amplitude
[M.max_dPdt,ip] = max(dPdt); M.tPpos=t(ip);
[M.min_dPdt,in] = min(dPdt); M.tPneg=t(in);
[M.max_dZdt,ip] = max(dZdt); M.tZpos=t(ip);
[M.min_dZdt,in] = min(dZdt); M.tZneg=t(in);

% peak lag
lag=[]; 
for i=1:numel(tPkP)
    tz = tPkZ(find(tPkZ>=tPkP(i),1,'first'));
    if ~isempty(tz), lag(end+1)=tz - tPkP(i); 
    end 
end
M.lag_mean = mean(lag);  M.lag_median = median(lag);

% time below 1/100 km^2
thr=0.01; dtloc=t(2)-t(1);
M.time_below_P = dtloc*sum(P<thr);
M.time_below_Z = dtloc*sum(Z<thr);
end

%%  PLOTS for problem 7
function fig = timeseries(t, P, Z, caseName)
fig = figure('Color','w'); hold on; grid on;

plot(t, P, 'Color',[0.1 0.45 0.85], 'LineWidth',1.8, 'DisplayName','Hare P(t)');
plot(t, Z, 'Color',[0.85 0.15 0.2],  'LineWidth',1.8, 'DisplayName','Lynx Z(t)');

xlabel('Time (years)','FontWeight','bold');
ylabel('Population density (/km^2)','FontWeight','bold');
title(sprintf('Time Series — %s', caseName));

hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end

% adaptive range
ymin = min([P Z]); ymax = max([P Z]); rangeY = ymax - ymin;
ylim([max(0,ymin-0.1*rangeY) ymax+0.1*rangeY]); xlim([0 t(end)]);

ax = gca; ax.FontSize = 12; ax.LineWidth = 1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'XColor','k','YColor','k','Box','off'); 
draw_black_border(ax); pbaspect([16 9 1]);
end

function fig = plot_phase_plane(P, Z, alpha, beta, gamma, epsi, P_eq, Z_eq, caseName)
fig = figure('Color','w'); hold on; grid on;

% trajectory
plot(P, Z, 'k-', 'LineWidth',1.8, 'DisplayName','Trajectory');

% nullclines
P_line = linspace(0, max(P)*1.1, 500);
colY = [0 0.7 0];   
colX = [1 0 1];    
Z_null = (alpha / gamma) * ones(size(P_line));
P_null = (beta  / (epsi * gamma)) * ones(size(P_line));
plot(P_line, Z_null, 'Color',colY, 'LineWidth',2, 'DisplayName','dP/dt=0');
plot(P_null, P_line, 'Color',colX, 'LineWidth',2, 'DisplayName','dZ/dt=0');

% stationary points
plot(0, 0, 'kp', 'MarkerSize',12, 'MarkerFaceColor','k', 'DisplayName','Stationary (0,0)');
plot(P_eq, Z_eq, 'ro', 'MarkerSize',8, 'MarkerFaceColor','r', 'DisplayName','Stationary (P^*,Z^*)');

% others
xl = xlabel('Prey P (hares/km^2)','FontWeight','bold');   set(xl,'Color','k');
yl = ylabel('Predator Z (lynx/km^2)','FontWeight','bold'); set(yl,'Color','k');
title(sprintf('Phase Plane — %s', caseName));
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end
xlim([0 max(P)*1.05]); ylim([0 max(Z)*1.05]);

% style
ax = gca; ax.FontSize = 12; ax.LineWidth = 1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'Box','off'); 
ax.XColor = colX;   
ax.YColor = colY;   
draw_black_border(ax); pbaspect([16 9 1]);
end

%% PLOTS for Problem 8
function fig = timeseries_log(t, P_lv, Z_lv, P_log, Z_log, ttl)
fig = figure('Color','w'); hold on; grid on;

plot(t, P_lv,  'LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (LV)');
plot(t, Z_lv,  'LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (LV)');
plot(t, P_log, '--','LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (Logistic)');
plot(t, Z_log, '--','LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (Logistic)');

xlabel('Time (years)','FontWeight','bold');     
ylabel('Population density (/km^2)','FontWeight','bold');
title(ttl);

hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end

% adaptive range
ymin = min([P_lv Z_lv P_log Z_log]); ymax = max([P_lv Z_lv P_log Z_log]); dy = ymax - ymin;
ylim([max(0,ymin-0.1*dy) ymax+0.1*dy]); xlim([0 t(end)]);

ax = gca; ax.FontSize = 12; ax.LineWidth = 1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'XColor','k','YColor','k','Box','off'); 
draw_black_border(ax); pbaspect([16 9 1]);
end

function fig = plot_phase_plane_logistic(P, Z, alpha, beta, gamma, epsi, K, P_eq, Z_eq, caseName)
fig = figure('Color','w'); hold on; grid on;

% trajectory
plot(P, Z, 'k-', 'LineWidth',1.8, 'DisplayName','Trajectory');

% nullclines 
u = linspace(0, max(max(P)*1.1, K), 500);
colY = [0 0.7 0];  
colX = [1 0 1];     
v_prey = (alpha/gamma) * (1 - u./K);                 
plot(u, v_prey, 'Color',colY, 'LineWidth',2, 'DisplayName','dP/dt=0 (logistic)');
plot((beta/(epsi*gamma))*ones(size(u)), u, 'Color',colX, 'LineWidth',2, 'DisplayName','dZ/dt=0');

% stationary points
plot(0, 0, 'kp', 'MarkerSize',12, 'MarkerFaceColor','k', 'DisplayName','Stationary (0,0)');
plot(P_eq, Z_eq, 'ro', 'MarkerSize',8, 'MarkerFaceColor','r', 'DisplayName','Stationary (P^*,Z^*)');

% others
xl = xlabel('Prey P (hares/km^2)','FontWeight','bold');   set(xl,'Color','k');
yl = ylabel('Predator Z (lynx/km^2)','FontWeight','bold'); set(yl,'Color','k');
title(sprintf('Phase Plane — %s', caseName));
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end

xlim([0 max([P u])*1.05]); ylim([0 max(Z)*1.05]);

% style
ax = gca; ax.FontSize = 12; ax.LineWidth = 1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'Box','off'); 
ax.XColor = colX;   
ax.YColor = colY;   
draw_black_border(ax); pbaspect([16 9 1]);
end

%% PLOTS for Problem 9 
function fig = timeseries_holling(t, P_lv, Z_lv, P_h, Z_h, ttl)
fig = figure('Color','w'); hold on; grid on;

plot(t, P_lv,  'LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (LV)');
plot(t, Z_lv,  'LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (LV)');
plot(t, P_h,  '--','LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (Holling-II)');
plot(t, Z_h,  '--','LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (Holling-II)');
xlabel('Time (years)','FontWeight','bold');
ylabel('Population density (/km^2)','FontWeight','bold');
title(ttl);
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end
ymin=min([P_lv Z_lv P_h Z_h]); ymax=max([P_lv Z_lv P_h Z_h]); dy=ymax-ymin;
ylim([max(0,ymin-0.1*dy) ymax+0.1*dy]); xlim([0 t(end)]);
ax=gca; ax.FontSize=12; ax.LineWidth=1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'XColor','k','YColor','k','Box','off'); draw_black_border(ax); pbaspect([16 9 1]);
end

function fig = phaseplane_holling(P, Z, alpha, beta, gamma, epsi, kappa, P_eq, Z_eq, ttl)
fig = figure('Color','w'); hold on; grid on;

% trajectory
plot(P, Z, 'k-', 'LineWidth',1.8, 'DisplayName','Trajectory');

% nullclines
u = linspace(0, max(P)*1.1, 500);
colY = [0 0.7 0];  
colX = [1 0 1];    
v_prey = (alpha/gamma) * (1 + gamma*kappa*u);
plot(u, v_prey, 'Color',colY, 'LineWidth',2, 'DisplayName','dU/dt=0 (Holling-II)');
if (epsi - beta*kappa) > 0
    Uv = beta / (gamma*(epsi - beta*kappa));
    plot(Uv*ones(size(u)), u, 'Color',colX, 'LineWidth',2, 'DisplayName','dV/dt=0');
else
    text(0.02*max(P), 0.9*max(Z), 'No positive dV/dt=0 nullcline', 'Color', colX);
end

% stationary points
plot(0, 0, 'kp', 'MarkerSize',12, 'MarkerFaceColor','k', 'DisplayName','(0,0)');
if ~isnan(P_eq) && ~isnan(Z_eq)
    plot(P_eq, Z_eq, 'ro', 'MarkerSize',8, 'MarkerFaceColor','r', 'DisplayName','(U^*,V^*)');
end

% others
xl = xlabel('Prey U (hares/km^2)','FontWeight','bold');   set(xl,'Color','k');
yl = ylabel('Predator V (lynx/km^2)','FontWeight','bold'); set(yl,'Color','k');
title(ttl);
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end
xlim([0 max(P)*1.05]); ylim([0 max(Z)*1.05]);

% style
ax = gca; ax.FontSize=12; ax.LineWidth=1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'Box','off'); ax.XColor=colX; ax.YColor=colY; draw_black_border(ax); pbaspect([16 9 1]);
end

%% PLOTS for Problem 10
function [P, Z] = simulate_logistic_holling(P0, Z0, alpha, beta, gamma, epsi, K, kappa, t, dt)
nT = numel(t); P=zeros(1,nT); Z=zeros(1,nT); P(1)=P0; Z(1)=Z0;
for k=1:nT-1
    func = (gamma*P(k)) / (1 + gamma*kappa*P(k)); 
    dP = alpha*P(k)*(1 - P(k)/K) - func*Z(k);
    dZ = epsi*func*Z(k) - beta*Z(k);
    P(k+1) = max(P(k) + dt*dP, 0);
    Z(k+1) = max(Z(k) + dt*dZ, 0);
end
end

function M = compute_M_logH(P, Z, t, alpha, beta, gamma, epsi, K, kappa)
func = (gamma.*P) ./ (1 + gamma*kappa.*P);
dPdt = alpha.*P.*(1 - P./K) - func.*Z;
dZdt = epsi.*func.*Z - beta.*Z;

M.P_min = min(P); M.P_max = max(P);
M.Z_min = min(Z); M.Z_max = max(Z);

burn_in = 5; mask = t>=burn_in;
[tPkP,~,iPkP] = detect_peaks(P(mask), t(mask), 1.0);
[tPkZ,~,iPkZ] = detect_peaks(Z(mask), t(mask), 1.0);
M.PerP = mean(diff(tPkP)); 
if isempty(M.PerP), M.PerP=NaN; end
M.PerZ = mean(diff(tPkZ)); 
if isempty(M.PerZ), M.PerZ=NaN; end
M.AmpP = mean_peak(P(mask), iPkP);
M.AmpZ = mean_peak(Z(mask), iPkZ);

[M.max_dPdt,ip] = max(dPdt); M.tPpos=t(ip);
[M.min_dPdt,in] = min(dPdt); M.tPneg=t(in);
[M.max_dZdt,ip] = max(dZdt); M.tZpos=t(ip);
[M.min_dZdt,in] = min(dZdt); M.tZneg=t(in);

% peak lag
lag=[]; 
for i=1:numel(tPkP)
    tz = tPkZ(find(tPkZ>=tPkP(i),1,'first'));
    if ~isempty(tz), lag(end+1)=tz - tPkP(i); 
    end
end
M.lag_mean = mean(lag);  
M.lag_median = median(lag);

thr=0.01; dtloc=t(2)-t(1);
M.time_below_P = dtloc*sum(P<thr);
M.time_below_Z = dtloc*sum(Z<thr);
end

function fig = timeseries_logH(t, P_lv, Z_lv, P_lh, Z_lh, ttl)
fig = figure('Color','w'); hold on; grid on;
plot(t, P_lv,  'LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (LV)');
plot(t, Z_lv,  'LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (LV)');
plot(t, P_lh, '--','LineWidth',1.8,'Color',[0.1 0.45 0.85],'DisplayName','Prey (Log+Hol)');
plot(t, Z_lh, '--','LineWidth',1.8,'Color',[0.85 0.15 0.2],'DisplayName','Pred (Log+Hol)');
xlabel('Time (years)','FontWeight','bold');
ylabel('Population density (/km^2)','FontWeight','bold');
title(ttl);
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end
ymin=min([P_lv Z_lv P_lh Z_lh]); ymax=max([P_lv Z_lv P_lh Z_lh]); dy=ymax-ymin;
ylim([max(0,ymin-0.1*dy) ymax+0.1*dy]); xlim([0 t(end)]);
ax=gca; ax.FontSize=12; ax.LineWidth=1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'XColor','k','YColor','k','Box','off'); draw_black_border(ax); pbaspect([16 9 1]);
end

function fig = plot_phase_plane_logH(P, Z, alpha, beta, gamma, epsi, K, kappa, Ueq, Veq, ttl)
fig = figure('Color','w'); hold on; grid on;

% trajectory
plot(P, Z, 'k-', 'LineWidth',1.8, 'DisplayName','Trajectory');

% nullclines
u = linspace(0, max(max(P)*1.1, K), 600);
colY = [0 0.7 0];  
colX = [1 0 1];    
v_u0 = (alpha/gamma) .* (1 + gamma*kappa.*u) .* (1 - u./K);
plot(u, v_u0, 'Color',colY, 'LineWidth',2, 'DisplayName','dU/dt=0 (Log+Hol)');
if (epsi - beta*kappa) > 0
    Uv = beta / (gamma*(epsi - beta*kappa));
    plot(Uv*ones(size(u)), u, 'Color',colX, 'LineWidth',2, 'DisplayName','dV/dt=0');
else
    text(0.02*max(P), 0.9*max(Z), 'No positive dV/dt=0 nullcline', 'Color',colX);
end

% stationary points
plot(0, 0, 'kp', 'MarkerSize',12, 'MarkerFaceColor','k', 'DisplayName','(0,0)');
if ~isnan(Ueq) && ~isnan(Veq)
    plot(Ueq, Veq, 'ro', 'MarkerSize',8, 'MarkerFaceColor','r', 'DisplayName','(U^*,V^*)');
end

% others
xl = xlabel('Prey U (hares/km^2)','FontWeight','bold');   set(xl,'Color','k');
yl = ylabel('Predator V (lynx/km^2)','FontWeight','bold'); set(yl,'Color','k');
title(ttl);
hide_name(gca, {'data1','data2','data3','data4'});
lgd = legend('Location','best'); if ~isempty(lgd)&&isvalid(lgd), lgd.AutoUpdate='off'; end
xlim([0 max([P u])*1.05]); ylim([0 max(Z)*1.05]);

% style
ax = gca; ax.FontSize=12; ax.LineWidth=1.2;
set(ax,'GridColor','k','MinorGridColor','k','XMinorGrid','on','YMinorGrid','on');
set(ax,'Box','off'); ax.XColor=colX; ax.YColor=colY; draw_black_border(ax); pbaspect([16 9 1]);
end

%% utilities 
function [tpk, ypk, ipk] = detect_peaks(y, t, min_sep_years)
% peak detection
if numel(y) < 3, tpk=[]; ypk=[]; ipk=[]; return; end
dt = t(2)-t(1);
cand = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
if isempty(cand), tpk=[]; ypk=[]; ipk=[]; return; end
min_sep = max(1, round(min_sep_years/dt));
sel = cand(1);
for i = 2:numel(cand)
    if cand(i) - sel(end) >= min_sep
        sel(end+1) = cand(i); 
    elseif y(cand(i)) > y(sel(end))
        sel(end) = cand(i);
    end
end
ipk = sel(:);
tpk = t(ipk);
ypk = y(ipk);
end

function amp = mean_peak(y, ipk)
% mean amplitude
if numel(ipk) < 2, amp = NaN; return; end
amps = zeros(numel(ipk)-1,1);
for i = 1:numel(ipk)-1
    seg = y(ipk(i):ipk(i+1));
    amps(i) = y(ipk(i)) - min(seg);
end
amp = mean(amps);
end

function draw_black_border(ax)
xl = xlim(ax); yl = ylim(ax);
hold(ax,'on');
line([xl(1) xl(2)], [yl(1) yl(1)], 'Color','k','LineWidth',1,'HandleVisibility','off'); % bottom
line([xl(1) xl(2)], [yl(2) yl(2)], 'Color','k','LineWidth',1,'HandleVisibility','off'); % top
line([xl(1) xl(1)], [yl(1) yl(2)], 'Color','k','LineWidth',1,'HandleVisibility','off'); % left
line([xl(2) xl(2)], [yl(1) yl(2)], 'Color','k','LineWidth',1,'HandleVisibility','off'); % right
set(ax,'Layer','top');
end

function hide_name(ax, names2hide)
for i = 1:numel(names2hide)
    h = findobj(ax,'-property','DisplayName','DisplayName',names2hide{i});
    for k = 1:numel(h)
        h(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
end

function print_delta(label, base, b, c)
fprintf('  • %-20s :  A=%.3f | B=%.3f (ΔB-A=%+.3f) | C=%.3f (ΔC-A=%+.3f)\n', ...
    label, base, b, b-base, c, c-base);
end

function print_cmp(label, a, b)
if numel(a)==2 && numel(b)==2
    fprintf('  %-22s : LV=[%g,%g] | Log=[%g,%g]\n', label, a(1), a(2), b(1), b(2));
else
    fprintf('  %-22s : LV=%g | Log=%g (Δ=%+g)\n', label, a, b, b-a);
end
end

function print_name(label, a, b, nameB)
if numel(a)==2 && numel(b)==2
    fprintf('  %-22s : LV=[%g,%g] | %s=[%g,%g]\n', label, a(1), a(2), nameB, b(1), b(2));
else
    fprintf('  %-22s : LV=%g | %s=%g (Δ=%+g)\n', label, a, nameB, b, b-a);
end
end
