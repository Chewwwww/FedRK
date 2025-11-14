% This code generates the results in section 3.1 of the paper (comparison of FedRK, FedAvg and FedProx). 

%% Set up the parameters of the experiments.
% m: number of rows
% n: number of columns
% T: number of federated rounds
% M: number of clients
% S: number of clients selected in each federated round
% tau: number of local iterations 

clear;

m = 2048;

n = 1024;

T = 1500; 

M = 16;
%M = 512;
%M = 1024

S = 5;


seed = RandStream('mlfg6331_64'); % random seed used in selection of subgroup of clients

%% Set up the system Ax = b

A = randn(m, n);

xstar = ones(n, 1);

b = A*xstar;

% generate a partition of the data 

p = partition(m, M, true); 

%% Run federated Kaczmarz (Algorithm 1). 

tau = [1, 50];

N = length(tau);

err_rk = zeros(T, N); % the distance between x and the sparse solution

for iter = 1: N

  %x = randn(n, 1); % the initial position
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

%% Run FedAvg 

tau = [1, 50];

N = length(tau);

err_avg = zeros(T, N); % the distance between x and the sparse solution

for iter = 1: N
  x = zeros(n, 1); % initialize the model

  err_avg(1, iter) = norm(x - xstar);

  xloc = zeros(n, S); % the positions of local clients

  for t = 1: T-1
    
    index = datasample(seed, 1:M, S, 'Replace', false);

    for i = 1: S
    
      xloc(:, i) = gd(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), tau(iter), 0.0005);
  
    end

    x = xloc*(p(index, 2) - p(index, 1)) / sum(p(index, 2) - p(index, 1));

    err_avg(t+1, iter) = norm(x - xstar) / norm(xstar);
  
  end
end

%% Run FedProx 

tau = [1, 50];

N = length(tau);

err_prox = zeros(T, N); % the distance between x and the sparse solution

for iter = 1: N
  x = zeros(n, 1); % initialize the model

  err_prox(1, iter) = norm(x - xstar);

  xloc = zeros(n, S); % the positions of local clients

  for t = 1: T-1
    
    index = datasample(seed, 1:M, S, 'Replace', false);

    for i = 1: S
    
      %xloc(:, i) = gd(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), tau(iter), 0.0005);
      xloc(:, i) = pgd(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), tau(iter), 0.0005, 0.25);
    end

    x = xloc*(p(index, 2) - p(index, 1)) / sum(p(index, 2) - p(index, 1));

    err_prox(t+1, iter) = norm(x - xstar) / norm(xstar);
  
  end
end

%% Plot the results

% Create figure
figure('Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');

hold on;

% Plot

lstyle = ["-", "--", ":", "-."];

for i = 1: N
  
  plot(log10(err_rk(:, i)), 'LineWidth', 2, 'LineStyle', lstyle(i));

end

for i = N: N
  
  plot(log10(err_avg(:, i)), 'LineWidth', 2, 'LineStyle', lstyle(3));

end

for i = N: N
  
  plot(log10(err_prox(:, i)), 'LineWidth', 2, 'LineStyle', lstyle(3));

end

% Labels and Title (with LaTeX)
xlabel('federated rounds', 'Interpreter', 'latex', 'FontSize', 18);

ylabel('$\log_{10}(\Vert x - x^* \Vert)$', 'Interpreter', 'latex', 'FontSize', 18);

title('Comparison of federated algorithms', 'Interpreter', 'latex', 'FontSize', 20);

% Legend
legend({'FedRK$(\tau=1)$', 'FedRK$(\tau=50)$', 'FedAvg$(\tau=50, \eta=0.0005)$', 'FedProx$(\tau=50, \eta=0.0005, \mu=0.25)$'}, ...
       'Interpreter', 'latex', ...
       'FontSize', 12, ...
       'Location', 'southwest');

% Axes settings
%grid on;
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.Box = 'on';

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save as high-resolution image
print('experiments_pic1','-dpng','-r600'); % Save as PNG, 600 dpi

