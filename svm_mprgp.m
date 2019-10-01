clear ALL;

addpath('mprgp')
addpath('problems/svm')

mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
rtol = 1e-1;

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
    
    % absolute tolerance
    atol = rtol*norm(b);
    
    % Lower bound
    mprgpctx.lb = zeros(size(y, 1), 1);
    
    %l1-loss
    fprintf('l1-loss\n')
    for k = 1:length(C)
        fprintf('C = %.f\n', C{k})
        mprgpctx.ub = C{k} * ones(size(y, 1), 1);

        % training
        tic;
        [x, flg, k] = mprgp(matvec, b, x_init, atol, 10000, mprgpctx);
        elapsed = toc;
        fprintf('%d; %d; %d; %d; %d; %e\n', flg, k(1), k(2), k(3), k(4), elapsed);
        
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
        [x, flg, k] = mprgp(A_reg, b, x_init, atol, 10000, mprgpctx);
        elapsed = toc;
        fprintf('%d; %d; %d; %d; %d; %e\n', flg, k(1), k(2), k(3), k(4), elapsed);

        % reconstruction of normal vector
        w = X' * Y * x;
        
        % performance scores of model
        scores = svm_eval(Test.X, Test.y, w);
        fprintf('%.2f%%; %.2f%%; %.2f%%; %.2f%%\n', scores(1), scores(2), scores(3), scores(3));
    end
end
