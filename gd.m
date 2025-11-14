function y = gd(A, x, b, k, eta)
  % A: the matrix 
  % x: the inital 
  % b: the right hand side
  % k: number of local iterations
  % eta: step size

  for i = 1: k

    x = x - eta * 0.5 * A'*(A*x - b);

  end

  y = x;
end