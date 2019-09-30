function [flg] = solved(nrmg,tol)
  flg = 0; % did not converge
  if nrmg <= tol
    flg = 1; % solver convergence
  end
end

