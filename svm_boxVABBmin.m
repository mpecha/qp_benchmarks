clear ALL;

addpath('boxABBmin_solvers')
addpath('problems/svm')

tol    = 1e-1;     % tolerance for the stopping criterion
tau    = 0.5;      % switching parameter
factor = 1.1;      %for factor=1 --> standard rule with fixed parameter tau
Mrho   = 3;        % memory for BB2
optioniterati = 0; % 1  -->  save iterates
epsi = 1e-7;

datasets = ["mushrooms"; "phishing"];
C = {1, 5, 10, 50, 100};

for d = 1:size(datasets,1)
    fprintf(strcat(datasets(d,:), '\n'))
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
    b = ones(size(y, 1), 1);
    % Initial guess
    x_init = zeros(size(y,1),1);
    
    % Lower bound
    lb = zeros(size(y, 1), 1);
    
    %l1-loss
    fprintf('l1-loss\n')
    for k = 1:length(C)
        fprintf('C = %.f\n', C{k})
        ub = C{k} * ones(size(y, 1), 1);
        
        % training
        tic;
        [x,info,f,step,err,vet_ls,iterati,times,nproj,nhprod] = ...
            boxABB_min_quad(A,b,lb,ub,x_init,tol,epsi,...
            'MAXIT',10000,...
            'VERBOSE',0,...
            'M', 10,...
            'OPTION',optioniterati,...
            'SAVE',1,...
            'STOPCRITERION',2,...
            'PARAMETER',[tau,factor,Mrho]);
        elapsed = toc;
        iter = length(step);
        fprintf('%d  k=%d \n',nhprod,iter);
        fprintf('Elapsed time = %e\n',elapsed)
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%;\n', scores(1), scores(2), scores(3), scores(3));
    end
    
    % l2-loss
    for k = 1:length(C)
        mprgpctx.ub = Inf * ones(size(y, 1), 1);
        A_reg = A + 1 / C{k} * eye(size(y, 1));
        matvec = @(x) A_reg*x;
        
        % training
        tic;
        [x,info,f,step,err,vet_ls,iterati,times,nproj,nhprod] = ...
            boxABB_min_quad(A,b,lb,ub,x_init,tol,epsi,...
            'MAXIT',10000,...
            'VERBOSE',0,...
            'M', 10,...
            'OPTION',optioniterati,...
            'SAVE',1,...
            'STOPCRITERION',2,...
            'PARAMETER',[tau,factor,Mrho]);
        elapsed = toc;
        iter = length(step);
        fprintf('%d  k=%d \n',nhprod,iter);
        fprintf('Elapsed time = %e\n',elapsed)
        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%\n', scores(1), scores(2), scores(3), scores(3));
    end
end
