function [x] = mprgpProj(x,mprgpctx)
  x = min(max(x,mprgpctx.lb),mprgpctx.ub);
end
