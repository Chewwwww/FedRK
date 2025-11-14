function y = fedrk(x, dloc, tau, option, varargin)
  if option == 'FR'
    
    if length(varargin) ~= 0
      error('incorrect number of arguments');
    end 

    % server RK

    S = size(dloc, 2);

    xloc = dloc + x*ones(1, S);

    b = zeros(S, 1);

    for i = 1: S
      b(i) = dloc(:, i)'*xloc(:, i);
    end
  
    y = rk(dloc', x, b, tau);

  elseif option == 'PT'
    
    if length(varargin) ~= 1
      error('incorrect number of arguments');
    end
  
    s = varargin{1};  % the argument specifies the sparsity level 

    % server RK

    S = size(dloc, 2);

    xloc = dloc + x*ones(1, S);

    b = zeros(S, 1);

    for i = 1: S
      b(i) = dloc(:, i)'*xloc(:, i);
    end
  
    y = rk(dloc', x, b, tau);

    % hard projection

    [u v] = sort(abs(y), 'descend');

    w = zeros(1, length(y))';
    v = v(1:floor(s));
    w(v) = 1;
    y = y.*w;
  
  end

end