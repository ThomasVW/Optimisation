function [p] = PenalisationDeCauchy(sigma, a, b, x, y)
  % Initialisation
  p = 0;
  N = size(x,1);

  for i=1:N
    % Calcul du r actuel.
    ri = a * x(i) + b - y(i);
    
    % Calcul de p et ajout au dernier p.
    p = p + (1/2)*log(1 + (ri/sigma)^2);
  end
end


