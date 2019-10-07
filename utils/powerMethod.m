function [ lambda, k] = powerMethod( matvec,x,tol,maxit )
%powerMethod - Estimate the largest eigenvalue
% https://en.wikipedia.org/wiki/Power_iteration
  lambda = norm(x);
  x = x/lambda;
  for k = 1:maxit
    lambda0 = lambda;
    x = matvec(x);
    lambda = norm(x);
    if abs(lambda-lambda0) < tol
      break;
    end
    x = x/lambda;
  end
end
