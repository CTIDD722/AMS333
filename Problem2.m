%% Problem 1 Yeast growth
% author: Yunpeng Chu
% Date:   2025/9/17
% Note:   Please use version 2013a or above

clear; clc; close all;
%% (a)
area   = 350;   % area
unit_K = 2.5e8; % carrying capacity per unit area
K      = area * unit_K;   

% 0.5(female ratio)*0.8(reproductive rate)*300(eggs/female)*0.2(egg->adult survival)
R0             = 0.5 * 0.8 * 300 * 0.2;     
                       
fprintf('===== (a) =====\n\n');
fprintf('K = %.3e (individuals)\n', K);
fprintf('R0 = %.2f \n\n', R0);

%% (b) 
% initialize parameters
b       = 1;
a       = R0 / K;                    
T_cycle = 17;  
N0      = 100;         
nmax    = 10;          
N       = zeros(1, nmax+1);
N(1)    = N0;
hassell = @(x) (R0 .* x) ./ (1 + a.*x).^b;

for n = 1:nmax
    N(n+1) = hassell(N(n));
end

% Population vs years
cycles     = 0:nmax;                 
N_ext      = [N, N(end)];
years     = cycles * T_cycle;
years_ext = [years, years(end)+T_cycle];

figure('Name','Time series','Color','w');
stairs(years_ext, N_ext, 'LineWidth', 1.8); hold on; grid on;
set(gca,'YScale','log');
plot(years, N, 'o', 'MarkerFaceColor','k', 'MarkerSize',5);
yline(K, '--', 'K');
xlabel('Years'); ylabel('N');
title('Population vs years');

% —— N(n+1) vs N(n)  
x = linspace(0, 1.2*K, 1200);
f = hassell(x);

figure('Name','(b3) N_{n+1} vs N_n','Color','w'); hold on; grid on;
plot(x, f, 'LineWidth', 1.8);       % Hassell 
plot(x, x, '--', 'LineWidth', 1.2); % y = x
% cobweb
Nc = N(1);
for k = 1:nmax
    plot([Nc Nc], [Nc hassell(Nc)], 'r-');          
    plot([Nc hassell(Nc)], [hassell(Nc) hassell(Nc)], 'r-'); 
    Nc = hassell(Nc);
end
plot(N(1:end-1), N(2:end), 'ko-', 'MarkerFaceColor','k');     
xlim([0, 1.2*K]); ylim([0, 1.2*K]);
xlabel('N_n'); ylabel('N_{n+1}');
title('N_{n+1} vs N_n (Hassell: b=1, a=R0/K)');
legend('Hassell','y=x','cobweb','trajectory','Location','NorthWest');

%% (c)
T_cycle = 17;           
alpha   = 0.99;           
Keff  = K * (R0 - 1) / R0;
hassell_step = @(N) (R0 * N) / (1 + a * N);

% find minimal n 
N0_test = 100;
N = N0_test;
n_need = 0;
while N < alpha * Keff
    N = hassell_step(N);
    n_need = n_need + 1;
    if n_need > 1e6, error('n too large—check parameters.'); end
end
T_need = n_need * T_cycle;

fprintf('===== Reaching %.0f%%·Keff with N0=%d =====\n', alpha*100, N0_test);
fprintf('Required cycles: %d  (~ %d years)\n\n', n_need, T_need);

% find minimal N0
function Nn_val = simulate_n_steps(N0, n, stepfun)
    Nn_val = N0;
    for ii = 1:n
        Nn_val = stepfun(Nn_val);
    end
end

% Bisection search 
function N0_min = find_min_N0(n, alpha, Keff, stepfun)
    tol = 1e-6;
    lo = 0;                
    hi = Keff;            
    while hi - lo > max(tol, 1e-12*max(1,hi))
        mid = 0.5*(lo + hi);
        if simulate_n_steps(mid, n, stepfun) >= alpha * Keff
            hi = mid;       
        else
            lo = mid;       
        end
    end
    N0_min = hi;
end

n_102 = 6; n_51 = 3;
N0_min_102y = find_min_N0(n_102, alpha, Keff, hassell_step);
N0_min_51y  = find_min_N0(n_51 , alpha, Keff, hassell_step);

fprintf('Within 102 years (6 cycles): N0 >= %.3e\n', N0_min_102y);
fprintf('Within  51 years (3 cycles): N0 >= %.3e\n\n', N0_min_51y);

%% (d) 
mass   = @(N) 0.002 * N;  
people = @(N) mass(N) / 70;

K_mass_kg = mass(K);
fprintf('===== (d) =====\n\n');
fprintf('K mass: %.3e kg  (~ %.1f people, assuming 70kg/person)\n', K_mass_kg, people(K));

N0_51_mass  = mass(N0_min_51y);
N0_102_mass = mass(N0_min_102y);
fprintf('51-year plan minimum release: %.2f tons (~ %.0f people)\n', N0_51_mass/1000, people(N0_min_51y));
fprintf('102-year plan minimum release: %.2f kg  (~ %.2f people)\n\n', N0_102_mass, people(N0_min_102y));
