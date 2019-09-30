function [gf,gc,gr] = mprgpSplit(x,g,mprgpctx)
  bounddiffl = x - mprgpctx.lb;
  bounddiffu = x - mprgpctx.ub;
  activesetl = (abs(bounddiffl) <= mprgpctx.settol);
  activesetu = (abs(bounddiffu) <= mprgpctx.settol);
  freeset = ~(activesetl+activesetu);
  gf = freeset.*g; % free
  gc = min(activesetl.*g,0.0); % chopped
  gc = gc + max(activesetu.*g,0.0); % chopped

  gset = (gf > 0.0);
  gr = min(gset.*bounddiffl/mprgpctx.abar,gf); % free reduced
  gr = gr + max(~gset.*bounddiffu/mprgpctx.abar,gf); % free reduced
end
