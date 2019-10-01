function [afeas] = mprgpFeas(x,p,mprgpctx)
      bounddiffl = x - mprgpctx.lb;
      bounddiffu = x - mprgpctx.ub;
      afeas = Inf;
      for i=1:length(x)
        %if (p(i) > mprgpctx.settol)
        if (p(i) > 0 && mprgpctx.lb(i) > -Inf)
          af = bounddiffl(i)/p(i);
          if af < afeas
            afeas = af;
          end
        elseif (p(i) < 0 && mprgpctx.ub(i) < Inf)
          af = bounddiffu(i)/p(i);
          if af < afeas
            afeas = af;
          end
        end
      end
end
