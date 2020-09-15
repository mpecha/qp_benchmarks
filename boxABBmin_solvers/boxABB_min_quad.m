function [x,info,f,a,err,vet_ls,iterati,times,nproj,nhprod] = ...
    boxABB_min_quad(A,b,lb,ub,x,tol,epsi,varargin)
% gradient projection method with steplength chosen in accordance
% with the BoxABB_min rule
%
%              min      1/2 *x'*A*x   -  b'x'
%              s. t.   lb <= x <= ub
%
% test for number of required parametres
if (nargin-length(varargin)) ~= 7
    error('Wrong number of required parameters');
end

% start the clock
t0 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%
% BoxABB_min default parameters
%%%%%%%%%%%%%%%%%%%%%%%%
maxit = 1000;          % maximum number of iterations
alphaini = 1.0;        % initial value for alpha
alphamin = 1e-10;      % alpha lower bound
alphamax = 1e6;		   % alpha upper bound
M = 10;                % memory in obj. function value (if M = 1 monotone)
gamma    = 1e-4;       % parameter for Armijo rule
beta     = 0.5;        % parameter for Armijo rule
Mrho     = 3;          % number of previous BB2 steplengths
tau      = 0.6;
factor   = 1;
%
option = 1;            % 1 ->  save iterates
verb = 0;              % 0 -> silent
stopcrit = 1;          % 1 -> norm(d,inf) <= tol
errcrit  = 0;          % 0 -> the exact solution is unknown
saveoption = 0;        % 0 -> save steplengths and function values
threshold = 10*eps;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SOL'
                sol = varargin{i+1};
                errcrit = 1;
            case 'MAXIT'
                maxit = varargin{i+1};
            case 'INITALPHA'
                alphaini = varargin{i+1};
            case 'ALPHAMIN'
                alphamin = varargin{i+1};
            case 'ALPHAMAX'
                alphamax = varargin{i+1};
            case 'GAMMA'
                gamma = varargin{i+1};
            case 'BETA'
                beta = varargin{i+1};
            case 'M'
                M = varargin{i+1};
            case 'STOPCRITERION'
                stopcrit = varargin{i+1};
            case 'OPTION'
                option = varargin{i+1};
            case 'VERBOSE'
                verb = varargin{i+1};
            case 'SAVE'
                saveoption =varargin{i+1};
            case 'PARAMETER'
                tau   = varargin{i+1}(1);
                factor= varargin{i+1}(2);
                Mrho  = varargin{i+1}(3);
            otherwise
                error('Unknown initialization option.');
        end
    end
end

nproj = 0;
nhprod = 0;

x = x(:);
x= min(max(x,lb),ub); nproj = nproj + 1;

if saveoption
    a = zeros(1,maxit);
    f = zeros(1,maxit+1);
else
    a = []; f = [];
end

if option
    iterati = cell(maxit+1,1);
    iterati{1} = x;
else
    iterati = [];
end

% error
if errcrit
    err = zeros(1, maxit+1);
    normsol = norm(sol);
    err(1) = norm(x-sol)/normsol;
else
    err = [];
end
times = zeros(1,maxit+1);


rho = alphaini;
g = A(x) - b; nhprod = nhprod + 1;
%%%% initial check of KKT
if stopcrit == 2
    phi0 = g;%.*(x~=lb & x~=ub) + min(0,g).*(x==lb)+ max(0,g).*(x==ub);
    normphi0 = norm(phi0);
    
    if verb >0
        fprintf('initial KKT: norm(phi0)=%g\n',normphi0);
    end
end


gold = g;
Fold = -Inf*ones(M,1);

vetrho = ones(Mrho,1)*alphamax;

vet_ls = [];
fnew = 0.5*(g-b)'*x;

if saveoption
    f(1) = fnew;
end

%%%%%%%%%%%%%%
% main loop
%%%%%%%%%%%%%
loop = true;
i = 0;
while loop
    
    vetrho(1:Mrho-1) = vetrho(2:Mrho);
    %%% descent direction
    d = min(max(x- rho*g,lb),ub); nproj = nproj + 1;
    d = d - x;
    
    if norm(d,inf) < epsi
        info = 0;
        break
    end
    
    i = i + 1;
    Fold(1:M-1) = Fold(2:M);
    Fold(M) = fnew;
    FR = max(Fold);
    xold = x;
    
    
    % nonmonotone linesearch GLL
    alpha = 1;
    gf = g'*d;
    xnew = x + alpha*d;
    Ad = A(d);  nhprod = nhprod + 1;
    gnew = g + alpha*Ad;
    fnew = 0.5*(gnew-b)'*xnew;
    infosearch = 1;
    while alpha > eps*norm(d)
        
        if FR - fnew < -gamma * alpha  * gf
            % fprintf('\n gd = %24.16e FR+sigma*alpha*gf=%g fnew=%g\n',gf, FR+sigma*alpha*gf, feval(fun,x+alpha*d));
            %  pause
            alpha = alpha * beta;
            xnew = x + alpha*d;
            gnew = g + alpha*Ad;
            fnew = 0.5*(gnew - b)'*xnew;
        else
            infosearch=0;
            break
        end
        
    end
    
    
    if infosearch>0
        fprintf(' \n bad behaviour in armijo: alpha<=eps*s %g  %g\n ',alpha,eps*norm(d));
    end
    
    if alpha < 1
        vet_ls = [vet_ls i];
    end
    
    
    s = alpha * d;
    x = xnew;
    
    if saveoption
        a(i)  = rho;
        f(i+1) = fnew;
    end
    if option
        iterati{i+1}=x;
    end
    if errcrit
        err(i+1) = norm(x-sol)/normsol;
    end
    
    g = gnew;
    y = g - gold;
    sy = s'*y;
    
    %%% MODIFIED BB2 in accordance with BoxBB2 updating rule
    y = y.* ((xold>lb | x>lb) & (xold<ub | x<ub));
    
    
    gold = g;
    
    if sy<0
        BB1 = min(10*rho,alphamax);
        BB2 = BB1;
        fprintf('negative curvature: iteration=%d',i);
        
    else
        BB1 = (s'*s)/sy;
        BB2 =  sy/(y'*y);
        BB1 = max(alphamin, min(alphamax,BB1));
        BB2 = max(alphamin, min(alphamax,BB2));
    end
    
    vetrho(Mrho) = BB2;
    
    if (BB2/BB1 < tau)
        rho = min(vetrho);
        tau = tau/factor;
    else
        rho = BB1;
        tau = tau*factor;
    end
    
%     
    gf = g.*(~(abs(x-lb)<=threshold | abs(x-ub)<=threshold)); 
    phi = gf + min(0,g).*(x==lb)+ max(0,g).*(x==ub);
    normphi = norm(phi);
    if verb>1
        fprintf('%d) norm(d)=%g rho=%g normKKT=%g', i, norm(d,inf),rho,normphi );
        if errcrit
            fprintf(' err=%g ', err(i+1));
        end
        fprintf('\n');
    end
    
    %%% stop criteria
    switch stopcrit
        case 1
            loop = (norm(s,inf) > tol | normphi>0.1) & (i< maxit) ;
            if verb>0
                fprintf('it=%d norm(d,inf)=%g \n', i,norm(d,inf));
            end
        case 2
            loop = (normphi > tol*normphi0)& (i< maxit);
            if verb>0
                fprintf('it=%d norm(phi)/norm(phi0)=%g \n',i, norm(phi)/normphi0);
            end
    end
    
    info = 0;
    times(i) = toc(t0);
end

iter = i;
if i>= maxit
    info=1;
end
times = times(1:iter);
if saveoption
    f = f(1:iter+1);
    a = a(1:iter);
end
if errcrit
    err = err(1:iter+1);
end
