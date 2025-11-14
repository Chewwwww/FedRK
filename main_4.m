% This code generates the results in section 4.3 of the paper.

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

m = 67;

n = 9;

s = 5;

T = 2000; 

M = 7;

S = 3;

N = 1;

%% read the data and perform the standaradization of features

data = readtable('data.csv');

X = table2array(data(:, 2: 9));

y = table2array(data(:, 10));

X = normalize(X);

%% perform train-test split

idx_train = strcmpi(data.train, {'T'});

idx_test = strcmpi(data.train, {'F'});

X_train = X(idx_train, :);

X_test = X(idx_test, :);

y_train = y(idx_train, :);

y_test = y(idx_test, :);


%% set up A b

A = [ones(size(X_train, 1), 1), X_train];

b = y_train;

%% Execute federated Kaczmarz with periodic thresholding at the server (Algorithm 2). N experiments are done in total.

orbit = zeros(9, T-1);

count = zeros(9, 1);

seed = RandStream('mlfg6331_64'); % random seed used in selection of subgroup of clients

for i = 1: N

% generate a partition of the data 

p = partition(m, M, true); 

x = randn(n, 1); % the initial position

xloc = zeros(n, S); % the positions of local clients

dloc = zeros(n, S); % the displacements of local clients

for t = 1: T-1
    
  index = datasample(seed, 1:M, S, 'Replace', false);

  for i = 1: S
    
    xloc(:, i) = rk(A(p(index(i), 1): p(index(i), 2), :), x, b(p(index(i), 1): p(index(i), 2), :), 5);

    dloc(:, i) = xloc(:, i) - x;
  
  end

  x = fedrk(x, dloc, 20, 'PT', s);

  orbit(:, t) = (x ~= 0);
  
end

end

count = sum(orbit, 2);

%% Plot the results

% Create figure
figure('Units', 'inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');

hold on;

% Plot
stem(count, "MarkerFaceColor", "red", "MarkerEdgeColor", "blue", "LineWidth", 2)

% Labels and Title (with LaTeX)
xlabel('features', 'Interpreter', 'latex', 'FontSize', 18);

ylabel('counts', 'Interpreter', 'latex', 'FontSize', 18);

title('feature selection', 'Interpreter', 'latex', 'FontSize', 20);

% Axes settings
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.Box = 'on';

% Tight layout
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Save as high-resolution image
print('experiments_pic4','-dpng','-r600'); % Save as PNG, 600 dpi
