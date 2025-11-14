% This code generates the results in section 3.1 of the paper (large system). 

%% Set up the parameters of the experiments.
% n: number of columns
% T: number of federated rounds
% M: number of clients
% S: number of clients selected in each federated round
% tau: number of local iterations 

clear;

n = 20480;

T = 1000; 

M = 10;

S = 5;

seed = RandStream('mlfg6331_64'); % random seed used in selection of subgroup of clients

%% Set up the system Ax = b

A = randn(2^5, n);

p = [1, 2^5];

for i = 1: M-1
    A = cat(1, A, (i+1)*randn(2^(5+i), n));

    p = cat(1, p, [p(i, 2)+1, p(i, 2)+2^(5+i)]);
end

xstar = ones(n, 1);

b = A*xstar;


%% Run federated Kaczmarz (Algorithm 1). 

tau = [50];

N = length(tau);

err_rk = zeros(T, N); % the distance between x and the sparse solution

for iter = 1: N

  x = zeros(n, 1);

  err_rk(1, iter) = norm(x - xstar);

  xloc = zeros(n, S); % the positions of local clients

  dloc = zeros(n, S); % the displacements of local clients

  for t = 1: T-1
    
    index = datasample(seed, 1:M, S, 'Replace', false);

    for i = 1: S
    
      xloc(:, i) = rk(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), tau(iter));

      dloc(:, i) = xloc(:, i) - x;
  
    end

    x = fedrk(x, dloc, 20, 'FR');

    err_rk(t+1, iter) = norm(x - xstar) / norm(xstar);
  
  end

end

%% Plot the results

% Create figure
figure('Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');

hold on;

% Plot

lstyle = ["-"];

for i = 1: N
  
  plot(log10(err_rk(:, i)), 'LineWidth', 2, 'LineStyle', lstyle(i));

end

% Labels and Title (with LaTeX)
xlabel('federated rounds (global iterations)', 'Interpreter', 'latex', 'FontSize', 18);

ylabel('$\log_{10}(\Vert x - x^* \Vert)$', 'Interpreter', 'latex', 'FontSize', 18);

title('FedRK for a large system', 'Interpreter', 'latex', 'FontSize', 20);

% Legend
legend({'FedRK$(\tau=50)$'}, ...
       'Interpreter', 'latex', ...
       'FontSize', 12, ...
       'Location', 'southwest');

% Axes settings

ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.Box = 'on';

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save as high-resolution image
print('experiments_pic1c','-dpng','-r600'); % Save as PNG, 600 dpi