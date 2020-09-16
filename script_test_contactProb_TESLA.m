clc;
clear all;
close all;

addpath('mprgp')
addpath('utils')
addpath('boxABBmin_solvers')
addpath('problems')

% MPRGP parameters
mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;


% GP-BoxVABBmin parameters
tau    = 0.5;     % switching parameter
factor = 1.1;      % for factor=1 --> standard rule with fixed parameter tau
Mrho   = 3;        % memory for BB2
M_Armijo = 10;     % 1 --> monotone linesearch
optioniterati = 0; % 1  -->  save iterates
stopcrit = 2;      % choose the stopping criterion
epsi = 1e-8;

% parameters shared by the methods
rtol = 1e-4;
maxit = 30000;

% Choose problem
files = ["ex1_100";"ex2_100";"jbearing2_50_50"];
%files = ["ex1_1000";"ex2_1000";"ex1_5000";
%          "ex2_100";"ex2_1000';'ex2_5000';
%          "jbearing2_50_50";"jbearing2_100_100";"jbearing2_200_50";"jbearing2_400_25"];

fileid = fopen('contactProb-output-TESLA.txt','a');

for i = 1:size(files,1)
    
    load(strcat('problems/',files(i,:),'.mat'))
    fprintf('\n'); fprintf(fileid,'\n');
    fprintf('File %s\n',files(i,:)); fprintf(fileid,'File %s\n',files(i,:));
    
    mprgpctx.lb = Problem.lb;
    if exist('Problem.ub','var')
        mprgpctx.ub = Problem.ub;
        ub =  Problem.ub; % for BoxVABBmin
    else
        mprgpctx.ub = ones(size(Problem.lb,1),1)*Inf;
        ub = ones(size(Problem.lb,1),1)*Inf; % for BoxVABBmin
    end
    atol = rtol*norm(Problem.b);
    matvec = @(x) Problem.A*x;
    
    %% MPRGP
    tic;
    [x1,flg,k] = mprgp(matvec,Problem.b,zeros(size(Problem.A,1),1),atol,maxit,mprgpctx);
    elapsed = toc;
    fprintf('Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
    fprintf('MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
        flg,k(1),k(2),k(3),k(4),elapsed);
    
    ij = find(abs(x1-Problem.lb)<=mprgpctx.settol | abs(x1-ub)<=mprgpctx.settol);
    szAct = length(ij) ;
    fprintf('size of final active set = %g\n', szAct)
    
    
    fprintf(fileid,'Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
    fprintf(fileid,'MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
        flg,k(1),k(2),k(3),k(4),elapsed);
    fprintf(fileid,'size of final active set = %g\n', szAct);
    
    %% MPRGPp
    tic;
    [x2,flg,k] = mprgp2(matvec,Problem.b,zeros(size(Problem.A,1),1),atol,maxit,mprgpctx);
    elapsed = toc;
    fprintf('Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
    fprintf('MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
        flg,k(1),k(2),k(3),k(4),elapsed);
    
    ij = find(abs(x2-Problem.lb)<=mprgpctx.settol | abs(x2-ub)<=mprgpctx.settol);
    szAct = length(ij) ;
    fprintf('size of final active set = %g\n', szAct)
    
    fprintf(fileid,'Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
    fprintf(fileid,'MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
        flg,k(1),k(2),k(3),k(4),elapsed);
    fprintf(fileid,'size of final active set = %g\n', szAct);
    
    %% GP-BoxVABBmin - tau = 0.5
    
    [x3,info,f,step,err,vet_ls,iterati,times,nproj,nhprod] = ...
        boxABB_min_quad(matvec,Problem.b,Problem.lb,ub,...
        zeros(size(Problem.A,1),1),rtol,epsi,...
        'MAXIT',maxit,...
        'VERBOSE',0,...
        'M', M_Armijo,...
        'OPTION',optioniterati,...
        'SAVE',1,...
        'STOPCRITERION',stopcrit,...
        'PARAMETER',[tau,factor,Mrho]);
    
    iter = length(step);
    
    fprintf('Alg:\t\t\t info\t   nHess\t  nProj\t  nBack  nIter\t  Time\n');
    fprintf('BoxVABBmin - tau=%g\t %d\t %d\t   %d\t   %d\t  %d\t %g\n',tau,info,...
        nhprod,nproj,length(vet_ls),iter,times(end));
    
    ij = find(abs(x3-Problem.lb)<=mprgpctx.settol | abs(x3-ub)<=mprgpctx.settol);
    szAct = length(ij) ;
    fprintf('size of final active set = %g\n', szAct);
    
    fprintf(fileid,'Alg:\t\t\t info\t   nHess\t  nProj\t  nBack  nIter\t  Time\n');
    fprintf(fileid,'BoxVABBmin  %d\t %d\t   %d\t   %d\t  %d\t %g\n',info,...
        nhprod,nproj,length(vet_ls),iter,times(end));
    fprintf(fileid,'size of final active set = %g\n', szAct);
    
    fprintf('\n'); fprintf(fileid,'\n');
    %compute errors
    e1 = norm(x1-x2,inf);
    fprintf('|x1-x2|_inf = %g\n',e1); fprintf(fileid,'|x1-x2|_inf = %g\n',e1);
    e2 = norm(x1-x3, inf);
    fprintf('|x1-x3|_inf = %g\n',e2); fprintf(fileid,'|x1-x3|_inf = %g\n',e2);
    e3 = norm(x2-x3,inf);
    fprintf('|x2-x3|_inf = %g\n',e3); fprintf(fileid,'|x2-x3|_inf = %g\n',e3);
end
fclose all;
rmpath('mprgp')
rmpath('utils')
rmpath('boxABBmin_solvers')
rmpath('problems')