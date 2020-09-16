clc;
clear all;
close all;

addpath('mprgp')
addpath('utils')
addpath('boxABBmin_solvers')
addpath('problems/svm')

% MPRGP parameters
mprgpctx.abarmult = 1.9;%1.95;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;

% BoxVABBmin parameters
tau    = 0.5;      % switching parameter
factor = 1.1;      % for factor=1 --> standard rule with fixed parameter tau
Mrho   = 3;        % memory for BB2
M_Armijo = 10;     % 1 --> monotone linesearch
optioniterati = 0; % 1  -->  save iterates
stopcrit = 2;      % choice of the stopping criterion
epsi = 1e-8;

% % parameters shared by the methods
rtol = 1e-1;      % tolerance for the stopping criterion
maxit = 10000;    % maximum number of iterations


fileid = fopen('svm-output-TESLA.txt','a');

datasets =  'mushrooms';%['phishing'];
C = {1, 5, 10, 50, 100};

for d = 1:size(datasets,1)
    fprintf(strcat(datasets(d,:), '\n'));
    fprintf(fileid,strcat(datasets(d,:), '\n'));
    dataset = strcat('problems/svm/', datasets(d,:), '/');
    
    % load training and test dataset
    load(strcat(dataset, 'training.mat'));
    load(strcat(dataset, 'test.mat'));
    
    X = Train.X;
    y = Train.y;
    Y = diag(y);
    
    % Hessian
    A = Y * X * X' * Y';
    matvec = @(x) A*x;
    % RHS
    b = ones(size(y, 2), 1);
    % Initial guess
    x_init = zeros(size(y,2),1);
    
    % absolute tolerance
    atol = rtol*norm(b);
    
    % Lower bound
    mprgpctx.lb = zeros(size(y, 1), 1);
    
    %% l1-loss
    fprintf('l1-loss\n');  fprintf(fileid,'l1-loss\n');
    for k = 1:length(C)
        fprintf('C = %.f\n', C{k}); fprintf(fileid,'C = %.f\n', C{k});
        mprgpctx.ub = C{k} * ones(size(y, 1), 1);
        
        fprintf('Alg; Solved;\tnHess;\tnCG;\t nExp;\tnProp;\tTime\n');
        fprintf(fileid,'Alg; Solved;\tnHess;\tnCG;\t nExp;\tnProp;\tTime\n');
        
        %% MPRGP
        % training
        tic;
        [x, flg, kvec] = mprgp(matvec, b, x_init, atol, maxit, mprgpctx);
        elapsed = toc;
        
        fprintf('MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        fprintf(fileid,'MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
        
        %% MPRGPp
        % training
        tic;
        [x, flg, kvec] = mprgp2(matvec, b, x_init, atol, maxit, mprgpctx);
        elapsed = toc;
        fprintf('MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        fprintf(fileid,'MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
        %% BoxVABBmin
        % training
        [x3,info,f3,step,err,vet_ls,iterati,times,nproj,nhprod] = ...
            boxABB_min_quad(matvec,b,mprgpctx.lb,mprgpctx.ub,x_init,rtol,epsi,...
            'MAXIT',maxit,...
            'VERBOSE',0,...
            'M', M_Armijo,...
            'OPTION',optioniterati,...
            'SAVE',1,...
            'STOPCRITERION',stopcrit,...
            'PARAMETER',[tau,factor,Mrho]);
        iter = length(step);
        
        fprintf('Alg:\t\t info  nHess\t   nProj\t  nBack\t  nIter\t  Time\n');
        fprintf(fileid,'Alg:\t info  nHess\t   nProj\t  nBack\t  nIter\t  Time\n');
        
        fprintf('BoxVABBmin - tau=%g  %d\t %d\t  %d\t   %d\t  %d\t  %g\n',...
            tau, info, nhprod, nproj, length(vet_ls),iter,times(end));
        fprintf(fileid,'BoxVABBmin  %d\t  %d\t   %d\t   %d\t  %d\t  %g\n',...
            info,nhprod,nproj,length(vet_ls),iter,times(end));
        
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
    end
    
    %% l2-loss
    fprintf('\n');  fprintf(fileid,'\n');
    fprintf('l2-loss\n');  fprintf(fileid,'l2-loss\n');
    for k = 1:length(C)
        fprintf('C = %.f\n', C{k}); fprintf(fileid,'C = %.f\n', C{k});
        mprgpctx.ub = Inf * ones(size(y, 1), 1);
        A_reg = A + 1 / C{k} * eye(size(y, 1));
        matvec = @(x) A_reg*x;
        
        %% MPRGP
        % training
        tic;
        [x, flg, kvec] = mprgp(matvec, b, x_init, atol, maxit, mprgpctx);
        elapsed = toc;
        
        fprintf('MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        fprintf(fileid,'MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
        
        %% MPRGPp
        % training
        tic;
        [x, flg, kvec] = mprgp2(matvec, b, x_init, atol, maxit, mprgpctx);
        elapsed = toc;
        fprintf('MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        fprintf(fileid,'MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg, kvec(1), kvec(2), kvec(3),kvec(4), elapsed);
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
        %% BoxVABBmin
        % training
        [x3,info,f3,step,err,vet_ls,iterati,times,nproj,nhprod] = ...
            boxABB_min_quad(matvec,b,mprgpctx.lb,mprgpctx.ub,x_init,rtol,epsi,...
            'MAXIT',maxit,...
            'VERBOSE',0,...
            'M', M_Armijo,...
            'OPTION',optioniterati,...
            'SAVE',1,...
            'STOPCRITERION',stopcrit,...
            'PARAMETER',[tau,factor,Mrho]);
        iter = length(step);
        
        fprintf('Alg:\t\t info  nHess\t   nProj\t  nBack\t  nIter\t  Time\n');
        fprintf(fileid,'Alg:\t info  nHess\t   nProj\t  nBack\t  nIter\t  Time\n');
        
        fprintf('BoxVABBmin - tau=%g  %d\t %d\t  %d\t   %d\t  %d\t  %g\n',...
            tau, info, nhprod, nproj, length(vet_ls),iter,times(end));
        fprintf(fileid,'BoxVABBmin  %d\t  %d\t   %d\t   %d\t  %d\t  %g\n',...
            info,nhprod,nproj,length(vet_ls),iter,times(end));
        
        
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2),...
            scores(3), scores(3));
        fprintf(fileid,'%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1),...
            scores(2),scores(3), scores(3));
        
    end
end
