clear; clc; close all;

%% Initialization
N   = 1e7;      % Total population
r0  = 3.0;      % Basic reproduction number r0
tau = 10;       % Mean infectious period tau (days)
I0  = 4;        % Initial number of infectious individuals
T   = 400;      % Number of simulated days
per = 1e5;      % Convert to "per 100,000 people"

%% Problem 1
[r0, tau, I0] = deal(r0, tau, I0);
[t, S, I, R, new] = sir_discrete(N, r0, tau, I0, T);

% 1) Peak current infections
[I_peak, x_I_peak] = max(I);
day      = t(x_I_peak);
act_case = I_peak * per;

% 2) Peak new daily infections
[new_peak, x_new_peak] = max(new);
day_peak  = t(x_new_peak) + 1;
new_cases = new_peak * per;

% 3) Final epidemic size
final_hat = R(end);
final_per = final_hat * 100;

% 4) Hospitalization needs
hosp_frac = 0.10;
beds  = hosp_frac * act_case;
admit = hosp_frac * new_cases;

% 5) Deaths
IFR   = 0.015;
tot_I = final_hat * N;
tot_D = IFR * tot_I;

fprintf('Baseline discrete SIR simulation\n');
beta  = r0 / tau;
gamma = 1 / tau;
fprintf('beta = %.4f, gamma = %.4f, r0 = %.2f, tau = %.1f days\n', ...
        beta, gamma, r0, tau);

fprintf('\n(1) Peak current infections:\n');
fprintf('    Day since first cases: %d days\n', day);
fprintf('    Fraction infected: %.4f\n', I_peak);
fprintf('    Active cases at peak: %.2f per 100,000 people\n', act_case);

fprintf('\n(2) Peak new daily infections:\n');
fprintf('    Day since first cases: %d days\n', day_peak);
fprintf('    New infections at peak: %.4f (fraction of population)\n', new_peak);
fprintf('    New daily cases at peak: %.2f per 100,000 people\n', new_cases);

fprintf('\n(3) Final epidemic size:\n');
fprintf('    Fraction of population infected: %.4f (%.2f%%)\n', final_hat, final_per);

fprintf('\n(4) Hospitalization needs:\n');
fprintf('    Peak beds needed: %.2f per 100,000 people\n', beds);
fprintf('    Peak daily admissions: %.2f per 100,000 people\n', admit);

fprintf('\n(5) Deaths:\n');
fprintf('    Total infected: %.0f people\n', tot_I);
fprintf('    Total deaths:  %.0f people\n', tot_D);

% Curves
figure;
plot(t, S, 'LineWidth', 2);
hold on;
plot(t, I, 'LineWidth', 2);
plot(t, R, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Baseline SIR: r_0 = 3, \tau = 10 days, I_0 = 4');
legend({'S(t)','I(t)','R(t)'}, 'Location', 'best');
grid on;

%% Q3 Case (i)
r0_3a  = 2 * r0;
tau_3a = tau;
I0_3a  = I0;

[t_3a, S_3a, I_3a, R_3a] = sir_discrete(N, r0_3a, tau_3a, I0_3a, T);

figure;
plot(t_3a, S_3a, 'LineWidth', 2);
hold on;
plot(t_3a, I_3a, 'LineWidth', 2);
plot(t_3a, R_3a, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q3 (i) (r_0 = 6, \tau = 10, I_0 = 4)');
legend({'S(t)','I(t)','R(t)'}, 'Location', 'best');
grid on;

%% Q3 Case (ii)
r0_3b  = r0;
tau_3b = tau / 2;
I0_3b  = I0;

[t_3b, S_3b, I_3b, R_3b] = sir_discrete(N, r0_3b, tau_3b, I0_3b, T);

figure;
plot(t_3b, S_3b, 'LineWidth', 2);
hold on;
plot(t_3b, I_3b, 'LineWidth', 2);
plot(t_3b, R_3b, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q3 (ii) (r_0 = 3, \tau = 5, I_0 = 4)');
legend({'S(t)','I(t)','R(t)'}, 'Location', 'best');
grid on;

%% Q3 Case (iii)
r0_3c  = r0;
tau_3c = tau;
I0_3c  = 400;

[t_3c, S_3c, I_3c, R_3c] = sir_discrete(N, r0_3c, tau_3c, I0_3c, T);

figure;
plot(t_3c, S_3c, 'LineWidth', 2);
hold on;
plot(t_3c, I_3c, 'LineWidth', 2);
plot(t_3c, R_3c, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q3 (iii) (r_0 = 3, \tau = 10, I_0 = 400)');
legend({'S(t)','I(t)','R(t)'}, 'Location', 'best');
grid on;

%% Q4 Case (i)
r0_4a  = 1.5;
tau_4a = tau;
I0_4a  = I0;

[t_4a, S_4a, I_4a, R_4a] = sir_discrete(N, r0_4a, tau_4a, I0_4a, T);
final_R_4a = R_4a(end);

figure;
plot(t_4a, I_4a, 'LineWidth', 2);
hold on;
plot(t_4a, R_4a, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q4 (i): r_0 = 1.5, \tau = 10, I_0 = 4');
legend({'I(t)','R(t)'}, 'Location', 'best');
grid on;

fprintf('Q4 (i): r0 = 1.5\n');
fprintf('    Final infected fraction R = %.4f (%.2f%%)\n', ...
        final_R_4a, final_R_4a*100);

%% Q4 Case (ii)
r0_4b  = 1.0;
tau_4b = tau;
I0_4b  = I0;

[t_4b, S_4b, I_4b, R_4b] = sir_discrete(N, r0_4b, tau_4b, I0_4b, T);
final_R_4b = R_4b(end);

figure;
plot(t_4b, I_4b, 'LineWidth', 2); hold on;
plot(t_4b, R_4b, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q4 (ii): r_0 = 1.0, \tau = 10, I_0 = 4');
legend({'I(t)','R(t)'}, 'Location', 'best');
grid on;

fprintf('Q4 (ii): r0 = 1.0\n');
fprintf('    Final infected fraction R = %.8f (%.6f%%)\n', ...
        final_R_4b, final_R_4b*100);

%% Q4 Case (iii)
r0_4c  = 0.5;
tau_4c = tau;
I0_4c  = I0;

[t_4c, S_4c, I_4c, R_4c] = sir_discrete(N, r0_4c, tau_4c, I0_4c, T);
final_R_4c = R_4c(end);

figure;
plot(t_4c, I_4c, 'LineWidth', 2); hold on;
plot(t_4c, R_4c, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title('Q4 (iii): r_0 = 0.5, \tau = 10, I_0 = 4');
legend({'I(t)','R(t)'}, 'Location', 'best');
grid on;

fprintf('Q4 (iii): r0 = 0.5\n');
fprintf('    Final infected fraction R = %.8f (%.6f%%)\n', ...
        final_R_4c, final_R_4c*100);

%% Q5
V_list = [0.25, 0.50, 0.75];

for idx = 1:length(V_list)
    V = V_list(idx);    % Prior immunity proportion in this case
    I0_hat = I0 / N;    % Initial infectious fraction
    R0_hat = V;         % Initial immune fraction
    S0_hat = 1 - I0_hat - R0_hat;

    [t_5, S_5, I_5, R_5] = sir_discrete_IC(r0, tau, S0_hat, I0_hat, R0_hat, T);

    R_tot   = R_5(end);
    ever_I  = R_tot - V;

    fprintf('Q5: prior immunity V = %.0f%%\n', V*100);
    fprintf('    Initial S0 = %.6f, I0 = %.6e, R0 = %.2f\n', S0_hat, I0_hat, R0_hat);
    fprintf('    Final immune fraction R_total = %.4f (%.2f%%)\n', R_tot, R_tot*100);
    fprintf('    Fraction infected in this epidemic (R_total - V) = %.4f (%.2f%%)\n', ...
            ever_I, ever_I*100);

    figure;
    plot(t_5, I_5, 'LineWidth', 2); hold on;
    plot(t_5, R_5, 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Fraction of population');
    title(sprintf('Q5: V = %.0f%%%% immune, r_0 = 3, \\tau = 10 days, I_0 = 4', V*100));
    legend({'I(t)','R(t)'}, 'Location', 'best');
    grid on;
end

%% Q6
Vc_r0 = 1 - 1 / r0;
fprintf('Q6: Herd immunity threshold for r0 = %.1f\n', r0);
fprintf('    V_c = 1 - 1/r0 = %.4f (%.2f%% of population)\n', ...
        Vc_r0, Vc_r0*100);

% Check R_eff for different V
V_list = [0.25 0.50 0.75];
for k = 1:length(V_list)
    v = V_list(k);
    Reff0 = r0 * (1 - v);
    fprintf('    For V = %.0f%%, R_eff,0= %.3f\n', v*100, Reff0);
end

%% Q7
r0_reduced = 1.5;
Vc_r0_1p5  = 1 - 1 / r0_reduced;

fprintf('Q7: Herd immunity threshold for r0 = %.1f\n', r0_reduced);
fprintf('    V_c = 1 - 1/r0 = %.4f (%.2f%% of population)\n', ...
        Vc_r0_1p5, Vc_r0_1p5*100);

I0_hat_7 = I0 / N;
R0_hat_7 = Vc_r0_1p5;
S0_hat_7 = 1 - I0_hat_7 - R0_hat_7;

[t_7, S_7, I_7, R_7] = sir_discrete_IC(r0_reduced, tau, ...
                                       S0_hat_7, I0_hat_7, R0_hat_7, T);

[I7_max, idx7] = max(I_7);
fprintf('    With r0 = %.1f and V = V_c, max I(t) = %.6e at t = %d days\n', ...
        r0_reduced, I7_max, t_7(idx7));
fprintf('    (Max I(t) is only slightly above the initial I_0 fraction, so no large epidemic.)\n');

figure;
plot(t_7, I_7, 'LineWidth', 2); hold on;
plot(t_7, R_7, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Fraction of population');
title(sprintf('Q7: r_0 = %.1f, V = V_c = %.2f%%%%, \\tau = 10, I_0 = 4', ...
              r0_reduced, Vc_r0_1p5*100));
legend({'I(t)','R(t)'}, 'Location','best');
grid on;


%% Discrete SIR (no prior immunity)
function [t, S_hat, I_hat, R_hat, new_inf_hat] = ...
         sir_discrete(N, r0, tau, I0, T_max)
    beta  = r0 / tau;
    gamma = 1 / tau;

    I0_hat = I0 / N;
    S0_hat = 1 - I0_hat;
    R0_hat = 0;

    t = 0:T_max;

    S_hat       = zeros(1, T_max+1);
    I_hat       = zeros(1, T_max+1);
    R_hat       = zeros(1, T_max+1);
    new_inf_hat = zeros(1, T_max);

    S_hat(1) = S0_hat;
    I_hat(1) = I0_hat;
    R_hat(1) = R0_hat;

    for k = 1:T_max
        new = beta * S_hat(k) * I_hat(k);
        new_inf_hat(k) = new;

        S_hat(k+1) = S_hat(k) - new;
        I_hat(k+1) = I_hat(k) + new - gamma * I_hat(k);
        R_hat(k+1) = R_hat(k) + gamma * I_hat(k);
    end
end

%% Discrete SIR with arbitrary initial S0, I0, R0 (fractions)
function [t, S_hat, I_hat, R_hat, new_inf_hat] = ...
         sir_discrete_IC(r0, tau, S0_hat, I0_hat, R0_hat, T_max)
    beta  = r0 / tau;
    gamma = 1 / tau;

    t = 0:T_max;

    S_hat       = zeros(1, T_max+1);
    I_hat       = zeros(1, T_max+1);
    R_hat       = zeros(1, T_max+1);
    new_inf_hat = zeros(1, T_max);

    S_hat(1) = S0_hat;
    I_hat(1) = I0_hat;
    R_hat(1) = R0_hat;

    for k = 1:T_max
        new = beta * S_hat(k) * I_hat(k);
        new_inf_hat(k) = new;

        S_hat(k+1) = S_hat(k) - new;
        I_hat(k+1) = I_hat(k) + new - gamma * I_hat(k);
        R_hat(k+1) = R_hat(k) + gamma * I_hat(k);
    end
end
