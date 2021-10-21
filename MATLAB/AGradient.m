function [derivee] = AGradient(sigma, a, b, x, y)
  derivee = 0;
  N = size(x,1);
  for i=1:N
    derivee = derivee + (((a*x(i)+b-y(i))*x(i))/(sigma^2 + ((a*x(i)+b-y(i))^2)));
  end
end
