function [afeas] = mprgpFeas(x,p,mprgpctx)
      alpha_f1 = inf;
      if ~isempty(mprgpctx.lb)
          i = find(p > 0);
          if min(size(i)) > 0
              alpha_f1 = min((x(i)-mprgpctx.lb(i))./p(i));
          end
      end

      alpha_f2 = inf;
      if ~isempty(mprgpctx.ub)
          j = find(p < 0);
          if min(size(j)) > 0
              alpha_f2 = min((x(j)-mprgpctx.ub(j))./p(j));
          end
      end

      afeas = min(alpha_f1,alpha_f2);
end
