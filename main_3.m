% This code generates the results in section 4.2 of the paper. 

%% Set up the parameters of the experiments.
% m: number of rows
% n: number of columns (n-s is the number of noise columns we add)
% s: the original dimension
% T: number of federated rounds
% M: number of clients
% S: number of clients selected in each federated round
% N: number of experiments

clear;

m = 2048;

n = [256, 512, 1024, 2048];

s = 256;

N = 1;

T = 400; 

M = 16;

S = 5;

%% setting up the system Ax  b, where x* is the least squares solution to the original system
AA = randn(m, s);

b = randn(m, 1);

xstar = lsqr(AA, b);

%% solve

err = zeros(T, length(n)); % the distance between x and the sparse solution

seed = RandStream('mlfg6331_64'); % random seed used in selection of subgroup of clients

for iter = 1: length(n)
  % generate the augmented matrix

  A = [AA, randn(m, n(iter)-s)];

  % generate a partition of the data 

  p = partition(m, M, true); 

  x = randn(n(iter), 1); % the initial position

  err(1, iter) = norm(x(1: s, 1) - xstar);

  xloc = zeros(n(iter), S); % the positions of local clients

  dloc = zeros(n(iter), S); % the displacements of local clients

  for t = 1: T-1
    
    index = datasample(seed, 1:M, S, 'Replace', false);

    for i = 1: S
    
      xloc(:, i) = rk(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), 20);

      dloc(:, i) = xloc(:, i) - x;
  
    end

    x = fedrk(x, dloc, 20, 'FR');

    err(t+1, iter) = norm(x(1: s, 1) - xstar);
  
  end

end

%% Plot the results

% Create figure
figure('Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');

hold on;

% Plot

lstyle = ["-", "--", ":", "-."];

for i = 1: length(n)
  
  plot(err(:, i), 'LineWidth', 2, 'LineStyle', lstyle(i));

end

% Labels and Title (with LaTeX)
xlabel('iterations', 'Interpreter', 'latex', 'FontSize', 18);

ylabel('error', 'Interpreter', 'latex', 'FontSize', 18);

title('FedRK with inconsistent system', 'Interpreter', 'latex', 'FontSize', 20);

% Legend
legend({'$n=256$', '$n=512$', '$n=1024$', '$n=2056$'}, ...
       'Interpreter', 'latex', ...
       'FontSize', 12, ...
       'Location', 'northeast');

% Axes settings
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.Box = 'on';

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save as high-resolution image
print('experiments_pic3','-dpng','-r600'); % Save as PNG, 600 dpi

