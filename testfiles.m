addpath('mprgp')
addpath('utils')

mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
rtol = 1e-12;
maxit = 30000;

files = ["ex1_1000";"ex2_100";"jbearing2_50_50"];
%files = ["ex1_100";"ex1_1000";"ex1_5000";
%          "ex2_100";"ex2_1000';'ex2_5000';
%          "jbearing2_50_50";"jbearing2_100_100";"jbearing2_200_50";"jbearing2_400_25"];


fprintf('File; Alg; Solved; nHess; nCG; nExp; nProp; Time\n');

for i = 1:size(files,1)
  load(strcat('problems/',files(i,:),'.mat'))

  mprgpctx.lb = Problem.lb;
  if exist('Problem.ub','var')
    mprgpctx.ub = Problem.ub;
  else
    mprgpctx.ub = ones(size(Problem.lb,1),1)*Inf;
  end
    %mprgpctx.lb = -ones(size(Problem.lb,1),1)*Inf;
  atol = rtol*norm(Problem.b);
  matvec = @(x) Problem.A*x;

  tic;
  [x,flg,k] = mprgp(matvec,Problem.b,zeros(size(Problem.A,1),1),atol,maxit,mprgpctx);
  elapsed = toc;
  if files(i,1:2) == "ex"
    pts = 0:1/(length(x)-1):1
    xa = @(z) 15.*z.^2-15*z;
    xexact = mprgpProj(xa(pts)',mprgpctx);
    relerr = norm(xexact - x)/norm(xexact)
    plot(1:length(x),x,1:length(x),xexact,1:length(x),mprgpctx.lb);
    pause
  else
   relerr = 0;
  end
  fprintf('%s; MPRGP; %d; %e; %d; %d; %d; %d; %e\n',files(i,:),flg,relerr,k(1),k(2),k(3),k(4),elapsed);
  tic;
  [x2,flg,k] = mprgp2(matvec,Problem.b,zeros(size(Problem.A,1),1),atol,maxit,mprgpctx);
  elapsed = toc;
  fprintf('%s; MPRGPp; %d; %e; %d; %d; %d; %d; %e\n',files(i,:),flg,relerr,k(1),k(2),k(3),k(4),elapsed);
end

% BQP problem
bqp_dirs = ["BQP-15000/"];
bqp_dir = "BQP/";
bqp_range = 1:36;
bqp_range = 4:4;

for i = 1:size(bqp_dirs,1)
  for j = bqp_range
    file = strcat('problems/',bqp_dir,bqp_dirs(i,:),'BQP',num2str(j),'.mat');
    load(file)
    mprgpctx.lb = l;
    mprgpctx.ub = u;
    matvec = @(x) MatVetProduct(d,P,x);
    x0 = min(max(x0vect,l),u); % TODO try also with x0 = 0?
    atol = rtol*norm(matvec(x0)-c);

    tic;
    [x,flg,k] = mprgp(matvec,c,x0,atol,maxit,mprgpctx);
    elapsed = toc;
    fprintf('%s; MPRGP; %d; %d; %d; %d; %d; %e\n',strcat('BQP',num2str(j)),flg,k(1),k(2),k(3),k(4),elapsed);
    tic;
    [x2,flg,k] = mprgp2(matvec,c,x0,atol,maxit,mprgpctx);
    elapsed = toc;
    fprintf('%s; MPRGPp; %d; %d; %d; %d; %d; %e\n',strcat('BQP',num2str(j)),flg,k(1),k(2),k(3),k(4),elapsed);
  end
end

