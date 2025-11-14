function y = rk(A, x, b, k)
  % A: the matrix 
  % x: the inital 
  % b: the right hand side
  % k: number of local iterations
  
  m = size(A, 1);

  p = zeros(m, 1);

  for i = 1: m
    p(i) = norm(A(i, :))^2;
  end

  p = p/sum(p);

  for i = 1: k
    ind = find(mnrnd(1, p));
    atemp = A(ind, :);
    btemp = b(ind);

    x = x + ((btemp - atemp*x) / norm(atemp)^2)*atemp';
  end

  y = x;
