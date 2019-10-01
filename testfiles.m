addpath('mprgp')

mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
rtol = 1e-6;

files = ["ex1_100";"ex2_100";"jbearing2_50_50"];
%files = ["ex1_100";"ex1_1000";"ex1_5000";
%          "ex2_100";"ex2_1000';'ex2_5000';
%          "jbearing2_50_50";"jbearing2_100_100";"jbearing2_200_50";"jbearing2_400_25"];
for i = 1:size(files,1)
  load(strcat('problems/',files(i,:),'.mat'))
  
  mprgpctx.lb = Problem.lb;
  if exist('Problem.ub','var')
    mprgpctx.ub = Problem.ub;
  else
    mprgpctx.ub = ones(size(Problem.lb,1),1)*Inf;
  end
  atol = rtol*norm(Problem.b);
  
  tic;
  [x,flg,k] = mprgp(Problem.A,Problem.b,zeros(size(Problem.A,1),1),atol,30000,mprgpctx);
  elapsed = toc;
  fprintf('%s; %d; %d; %d; %d; %d; %e\n',files(i,:),flg,k(1),k(2),k(3),k(4),elapsed);
end
