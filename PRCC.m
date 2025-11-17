tic; clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1);format long g;

% Set a random seed for reproducibility
rng(13579);

% Initial conditions
n0 = 1; m0 = 1; c0 = 1; z0 = [n0; m0; c0];

% Time points
t = linspace(0, 1000, 100000);

% Define the problem
problem.num_vars = 9;
problem.names = {'r_1', 'r_2', '\alpha_1', '\alpha_2', 'h_1', 'h_2', '\beta_1', '\beta_2', '\mu'};
problem.bounds = [1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 1, 10; 0, 1];
% r1 = params(1); r2 = params(2); alpha1 = params(3); alpha2 = params(4);
% h1 = params(5); h2 = params(6); beta1 = params(7); beta2 = params(8); mu = params(9);

% Generate samples
param_values = lhsdesign(2^8, problem.num_vars);
% % Scale samples to the bounds
for i = 1:problem.num_vars
    param_values(:, i) = param_values(:, i) * (problem.bounds(i, 2) - problem.bounds(i, 1)) + problem.bounds(i, 1);
end
% Run the model
Y = zeros(size(param_values, 1), length(t));
for i = 1:size(param_values, 1)
    params = param_values(i, :);
    r1 = params(1); r2 = params(2); alpha1 = params(3); alpha2 = params(4);
    h1 = params(5); h2 = params(6); beta1 = params(7); beta2 = params(8); mu = params(9);
    [~, solution] = ode45(@(t, z) GMM(t, z, r1, r2, alpha1, alpha2, h1, h2, beta1, beta2, mu), t, z0);
    Y(i, :) = solution(:, end);  % Collecting the GMM sub-population
end

% Calculate PRCC
prcc_values = zeros(1, problem.num_vars);
for i = 1:problem.num_vars
    prcc_values(i) = corr(param_values(:, i), Y(:, end), 'Type', 'Spearman');
end
disp('PRCC values:'); % disp(prcc_values);
% Convert parameter names to categorical
param_names = categorical(problem.names);

% Plotting the PRCC values
figure;
bar(param_names, prcc_values, 'FaceColor', 'flat');
xlabel('Parameters', 'FontSize', 12); ylabel('Sensitivity Index', 'FontSize', 12);
% title('PRCC Sensitivity Analysis');
toc

% The WT-GMM model without spatial variations terms is defined as follows
function dzdt = GMM(t, z, r1, r2, alpha1, alpha2, h1, h2, beta1, beta2, mu)
    n = z(1); m = z(2); c = z(3);
    dndt = r1 * n * (1 - n - alpha1 * m);
    dmdt = r2 * m * (1 - n - alpha2 * m);
    dcdt = h1 * n + h2 * m - (mu + beta1 * n + beta2 * m) * c;
    dzdt = [dndt; dmdt; dcdt];
end
