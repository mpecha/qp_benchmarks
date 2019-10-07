function [x,flg,k,iter] = mprgp(matvec,b,x,tol,maxit,mprgpctx)
% solve .5x'Ax -x'b s.t. x >= lb
% 
% tol = abstol
%
% mprgpctx
% % mprgpctx.lb = bound
% % mprgpctx.abarmult = (0,2], length of expansion step = abarmult*||A^{-1}||
% % mprgpctx.propConst = proportioning const.
% % mprgpctx.settol = 10*eps, splitting tolerance
% % mprgpctx.infeastol = 0, infeasibility tolerance

  k = [0 0 0 0]; %Hessian, CG, Exp, Prop
  iter = 0;
  if  ~exist('mrpgpctx.abar', 'var')
    %mprgpctx.abar = mprgpctx.abarmult/eigs(matvec,length(b),1);
    [mprgpctx.abar,iters] = powerMethod(matvec,ones(length(b),1),1e-4,50);
    mprgpctx.abar = mprgpctx.abarmult*mprgpctx.abar;
    k = k + [iters 0 0 0];
  end
  x = mprgpProj(x,mprgpctx); % project to feasible set
  g = matvec(x) - b; k = k + [1 0 0 0];
  [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
  p = gf;
  flg = solved(norm(gf+gc),tol);
  while ~flg && iter < maxit
    if gc'*gc <= mprgpctx.propConst^2 *gr'*gf % proportional
      Ap = matvec(p); k = k + [1 0 0 0];
      pAp = p'*Ap;
      acg = g'*p/pAp;
      afeas = mprgpFeas(x,p,mprgpctx);
      if acg <= afeas % cg
        step = 'c';
        x = x - acg*p;
        g = g - acg*Ap;
        [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
        bcg = gf'*Ap/pAp;
        p = gf - bcg*p;
        k = k + [0 1 0 0];
      else % expansion
        step = 'e';
        x = x - afeas*p;
        g = g - afeas*Ap;
        [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
        %x = x - mprgpctx.abar*gr;
        x = x - mprgpctx.abar*gf;
        x = mprgpProj(x,mprgpctx); % project to feasible set
        g = matvec(x) - b; k = k + [1 0 1 0];
        [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
        p = gf;
      end
    else % proportioning
      step = 'p';
      p = gc;
      Ap = matvec(p); k = k + [1 0 0 1];
      pAp = p'*Ap;
      acg = g'*p/pAp;
      x = x - acg*p;
      g = g - acg*Ap;
      [gf,gc,gr] = mprgpSplit(x,g,mprgpctx);
      p = gf;
    end
    %fprintf("%s ",step)
    flg = solved(norm(gf+gc),tol);
    iter = iter + 1;
  end
end

