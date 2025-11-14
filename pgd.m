function y = pgd(A, x, b, k, eta, mu)
  % A: the matrix 
  % x: the inital 
  % b: the right hand side
  % k: number of local iterations
  % eta: step size

  x0 = x;

  for i = 1: k

    x = x - eta * 0.5 * A'*(A*x - b) - mu * (x - x0);

  end

  y = x;
  
end