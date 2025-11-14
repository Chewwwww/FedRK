function p = partition(m, M, uni)
  % generate a partition of [M] = {1,...,M} into m nonempty parts
  % opt specifies the distribution of the cardinality of each parts

  % p(i, :) is a two dimensional vector specifying the start and the end of ith client's data, 
% e.g. p(i, :) = [3, 7] means rows 3 to 7 (inclusive) of the system are i-th client's data.

  p = zeros(M, 2);

  q = floor(linspace(1, m+1, M+1));

  for i = 1: M
    
    p(i, 1) = q(i);

    p(i, 2) = q(i+1) - 1;

  end

%  if uni

%    s = m/M; % the size of each parts

%    for i = 1: M
%        p(i, :) = [1 + (i-1)*s, i*s];
%    end

%  else
    
%    c = [ sort(datasample(1: m-1, M-1, 'Replace', false)), m] + 0.5;

%    iter = 0;

%    for i = 1: M
%        iter = iter + 1;

%        p(i, 1) = iter;

%        while iter + 1 < c(i)
%            iter = iter + 1
%        end

%        p(i, 2) = iter;
%    end
%  end
%end