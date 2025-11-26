tic; clear all; close all; clc
set(0, 'defaultaxesfontsize', 10, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1);
format long g;

rng(13579);

% Initial conditions
n0 = 1; n10 = 0; m0 = 1; m10 = 0; c0 = 1; c10 = 0; z0 = [n0; n10; m0; m10; c0; c10];

% Time points
t = linspace(0, 100, 1000000);

% Define problem
problem.num_vars = 18;
problem.names = {'r_1', 'r_2', '\alpha_1', '\alpha_2', '\sigma_1', '\sigma_2', '\chi_1', '\chi_2', 'D_1', 'D_2', 'D_3', 'h_1', 'h_2', '\beta_1', '\beta_2', '\mu', 'u0', 'V_g'};
problem.bounds = [0, 1; 0, 1; 0, 1; 0, 1; 1, 10; 1, 10; 0, 1; 0, 1; 0, 1; 0, 1; 0, 1; 1, 10; 1, 10; 1, 10; 1, 10; 0, 1; 0, 1; 0, 10];

% Generate samples
param_values = lhsdesign(2^8, problem.num_vars);
for i = 1:problem.num_vars
    param_values(:, i) = param_values(:, i) * (problem.bounds(i, 2) - problem.bounds(i, 1)) + problem.bounds(i, 1);
end

% Run model and collect outputs
Y = zeros(size(param_values, 1), 6); % n, m, c at final time
for i = 1:size(param_values, 1)
    params = param_values(i, :);
    [~, solution] = ode15s(@(t, z) GMM(t, z, params), t, z0);
    Y(i, :) = solution(end, :); % final values of n, m, c
end

% Compute PRCC for each variable
prcc_values = zeros(problem.num_vars, 6);
for j = 1:6
    for i = 1:problem.num_vars
        prcc_values(i, j) = corr(param_values(:, i), Y(:, j), 'Type', 'Spearman');
    end
end

disp('PRCC matrix (rows=parameters, cols=[n m c]):');
disp(prcc_values);

% Select only the parameters of interest
selected_params = {'r_1','r_2','\alpha_1','\alpha_2','\sigma_1','\sigma_2',...
                   '\chi_1','\chi_2','D_1','D_2','D_3','h_1','h_2',...
                   '\beta_1','\beta_2','\mu','u0'};

% Find their indices in problem.names
[~, idx] = ismember(selected_params, problem.names);

% figure('Position',[100 100 1400 800]); % make figure larger

% PRCC for n
figure
bar(prcc_values(idx,1));
xticks(1:length(selected_params));
xticklabels(selected_params);
xtickangle(45);
set(gca,'FontSize',10);
ylabel('PRCC (n)');
% title('PRCC for n');

% PRCC for m
figure
bar(prcc_values(idx,3));
xticks(1:length(selected_params)); 
xticklabels(selected_params);
xtickangle(45); set(gca,'FontSize',10);
ylabel('PRCC (m)');
% title('PRCC for m');

% PRCC for c
figure
bar(prcc_values(idx,5));
xticks(1:length(selected_params));
xticklabels(selected_params);
xtickangle(45); set(gca,'FontSize',10);
ylabel('PRCC (c)');
% title('PRCC for c');


% ODE system
function dzdt = GMM(t, z, params)
    r1 = params(1); r2 = params(2); alpha1 = params(3); alpha2 = params(4);
    sigma1 = params(5); sigma2 = params(6); chi1 = params(7); chi2 = params(8);
    D1 = params(9); D2 = params(10); D3 = params(11);
    h1 = params(12); h2 = params(13); beta1 = params(14); beta2 = params(15); mu = params(16);
    u0 = params(17); Vg = params(18);
    n = z(1); n1 = z(2); m = z(3); m1 = z(4); c = z(5); c1 = z(6);
    dndt = n1;
    dn1dt = (((u0 - Vg + sigma1*n1 + sigma2*m1)*c1 - h1*n + h2*m + c*(mu + beta1*n + beta2*m))*((-sigma2*n)*(-chi2*m) - (-chi1*n)*(D2 - sigma2*m)) + (-(-sigma2*n)*(-D3) + (-chi1*n)*(-sigma2*c))*((u0 - Vg + chi2*c1 + sigma1*n1 + sigma2*m1)*m1 - r2*m*(1 - m - alpha2*n)) + ((D2 - sigma2*m)*(-D3) - (-chi2*m)*(-sigma2*c))*((u0 - Vg + chi1*c1 + sigma1*n1 + sigma2*m1)*n1 - r1*n*(1 - n - alpha1*m)))/((D1 - sigma1*n)*(D2 - sigma2*m)*(-D3) - (D1 - sigma1*n)*(-chi2*m)*(-sigma2*c) - (-sigma2*n)*(-sigma1*m)*(-D3) + (-sigma2*n)*(-chi2*m)*(-sigma1*c) + (-chi1*n)*(-sigma1*m)*(-sigma2*c) - (-chi1*n)*(D2 - sigma2*m)*(-sigma1*c));
    dmdt = m1;
    dm1dt = -(((D1 - sigma1*n)*(-chi2*m) - (-chi1*n)*(-sigma1*m))*((u0 - Vg + sigma1*n1 + sigma2*m1)*c1 - h1*n + h2*m + c*(mu + beta1*n + beta2*m)) + (-(D1 - sigma1*n)*(-D3) + (-chi1*n)*(-sigma1*c))*((u0 - Vg + chi2*c1 + sigma1*n1 + sigma2*m1)*m1 - r2*m*(1 - m - alpha2*n)) + ((-sigma1*m)*(-D3) - (-chi2*m)*(-sigma1*c))*((u0 - Vg + chi1*c1 + sigma1*n1 + sigma2*m1)*n1 - r1*n*(1 - n - alpha1*m)))/((D1 - sigma1*n)*(D2 - sigma2*m)*(-D3) - (D1 - sigma1*n)*(-chi2*m)*(-sigma2*c) - (-sigma2*n)*(-sigma1*m)*(-D3) + (-sigma2*n)*(-chi2*m)*(-sigma1*c) + (-chi1*n)*(-sigma1*m)*(-sigma2*c) - (-chi1*n)*(D2 - sigma2*m)*(-sigma1*c));
    dcdt = c1;
    dc1dt = (((D1 - sigma1*n)*(D2 - sigma2*m) - (-sigma2*n)*(-sigma1*m))*((u0 - Vg + sigma1*n1 + sigma2*m1)*c1 - h1*n + h2*m + c*(mu + beta1*n + beta2*m)) + (-(D1 - sigma1*n)*(-sigma2*c) + (-sigma2*n)*(-sigma1*c))*((u0 - Vg + chi2*c1 + sigma1*n1 + sigma2*m1)*m1 - r2*m*(1 - m - alpha2*n)) + ((-sigma1*m)*(-sigma2*c) - (D2 - sigma2*m)*(-sigma1*c))*((u0 - Vg + chi1*c1 + sigma1*n1 + sigma2*m1)*n1 - r1*n*(1 - n - alpha1*m)))/((D1 - sigma1*n)*(D2 - sigma2*m)*(-D3) - (D1 - sigma1*n)*(-chi2*m)*(-sigma2*c) - (-sigma2*n)*(-sigma1*m)*(-D3) + (-sigma2*n)*(-chi2*m)*(-sigma1*c) + (-chi1*n)*(-sigma1*m)*(-sigma2*c) - (-chi1*n)*(D2 - sigma2*m)*(-sigma1*c));
    dzdt = [dndt; dn1dt; dmdt; dm1dt; dcdt; dc1dt];
end

