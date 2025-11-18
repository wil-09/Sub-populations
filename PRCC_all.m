tic; clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1);
format long g;

rng(13579);

% Initial conditions
n0 = 1; m0 = 1; c0 = 1; z0 = [n0; m0; c0];

% Time points
t = linspace(0, 1000, 100000);

% Define problem
problem.num_vars = 9;
problem.names = {'r_1', 'r_2', '\alpha_1', '\alpha_2', 'h_1', 'h_2', '\beta_1', '\beta_2', '\mu'};
% problem.bounds = [1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 0, 1];
problem.bounds = [1, 10; 1, 10; 0, 1; 0, 1; 1, 10; 1, 10; 1, 10; 1, 10; 0, 1];

% Generate samples
param_values = lhsdesign(2^8, problem.num_vars);
for i = 1:problem.num_vars
    param_values(:, i) = param_values(:, i) * (problem.bounds(i, 2) - problem.bounds(i, 1)) + problem.bounds(i, 1);
end

% Run model and collect outputs
Y = zeros(size(param_values, 1), 3); % n, m, c at final time
for i = 1:size(param_values, 1)
    params = param_values(i, :);
    [~, solution] = ode45(@(t, z) GMM(t, z, params), t, z0);
    Y(i, :) = solution(end, :); % final values of n, m, c
end

% Compute PRCC for each variable
prcc_values = zeros(problem.num_vars, 3);
for j = 1:3
    for i = 1:problem.num_vars
        prcc_values(i, j) = corr(param_values(:, i), Y(:, j), 'Type', 'Spearman');
    end
end

disp('PRCC matrix (rows=parameters, cols=[n m c]):');
disp(prcc_values);

% Plot grouped bar chart
figure;
bar(prcc_values);
set(gca, 'XTickLabel', problem.names);
legend({'n', 'm', 'c'}, 'Location', 'Best');
ylabel('PRCC');
title('PRCC for n, m, and c at final time point');
toc

% ODE system
function dzdt = GMM(t, z, params)
    r1 = params(1); r2 = params(2); alpha1 = params(3); alpha2 = params(4);
    h1 = params(5); h2 = params(6); beta1 = params(7); beta2 = params(8); mu = params(9);
    n = z(1); m = z(2); c = z(3);
    dndt = r1 * n * (1 - n - alpha1 * m);
    dmdt = r2 * m * (1 - n - alpha2 * m);
    dcdt = h1 * n + h2 * m - (mu + beta1 * n + beta2 * m) * c;
    dzdt = [dndt; dmdt; dcdt];
end