% This code generates the results in section 4.1 of the paper. 

%% Set up the parameters of the experiments.
% m: number of rows
% n: number of columns
% s: sparsity level
% sigma: standard deivation of the noise
% T: number of federated rounds
% M: number of clients
% S: number of clients selected in each federated round
% N: number of experiments

clear;

m = 256;

n = 1024;

s = 9;

sigma = 1e-2;

N = 10;

T = 1500; 

M = 16;

S = 5;

%% Set up the system Ax = Ax* + e, where x* is the target sparse solution and e is the noise.

A = randn(m, n);

xstar = rand(n, 1);

nonSupp = randsample(n, n-s);

xstar(nonSupp) = 0;

xstar = xstar/norm(xstar);

b = A*xstar + sigma*randn(m, 1);

%% Execute federated Kaczmarz with periodic thresholding at the server. N experiments are done in total.

count = zeros(n, 1);

err = zeros(T, N); % the distance between x and the sparse solution

seed = RandStream('mlfg6331_64'); % random seed used in selection of subgroup of clients

% generate a partition of the data 

p = partition(m, M, true); 

for iter = 1: N

  x = randn(n, 1); % the initial position

  err(1, iter) = norm(x - xstar);

  xloc = zeros(n, S); % the positions of local clients

  dloc = zeros(n, S); % the displacements of local clients

  for t = 1: T-1
    
    index = datasample(seed, 1:M, S, 'Replace', false);

    for i = 1: S
    
      xloc(:, i) = rk(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), 20);

      dloc(:, i) = xloc(:, i) - x;
  
    end

    x = fedrk(x, dloc, 20, 'PT', s);

    err(t+1, iter) = norm(x - xstar);
  
  end

  count = count + (x ~= 0);

end

%% Plot the results

% Create figure
figure('Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');

hold on;

% Plot
for i = 1: N
  
  plot(log10( err(:, i) ), 'LineWidth', 2);

end
%stem(count, "MarkerFaceColor", "red", "MarkerEdgeColor", "blue", "LineWidth", 2);

% Labels and Title (with LaTeX)
xlabel('federated rounds (global iterations)', 'Interpreter', 'latex', 'FontSize', 18);

ylabel('$\log_{10}(\Vert x - x^* \Vert)$', 'Interpreter', 'latex', 'FontSize', 18);

title('FedRK with hard thresholding at the server', 'Interpreter', 'latex', 'FontSize', 20);

% Axes settings
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.Box = 'on';

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save as high-resolution image
print('experiments_pic2a','-dpng','-r600'); % Save as PNG, 600 dpi

