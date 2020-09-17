clc;
clear ALL;
close all;

addpath('mprgp')
addpath('utils')
addpath('boxABBmin_solvers')

%  Create folders for saving outputs of each method
saveout = 1;
if saveout
    for j = 1:3
        folder = sprintf('output-data-%d',j);
        if ~exist(folder,'dir')
            mkdir(folder);
        end
    end
end

% MPRGP parameters
mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;

% BoxVABBmin parameters
tau    = 0.5;     % switching parameter
factor = 1.1;      % for factor=1 --> standard rule with fixed parameter tau
Mrho   = 3;        % memory for BB2
M_Armijo = 10;     % 1 --> monotone linesearch
optioniterati = 0; % 1  -->  save iterates
stopcrit = 2;      % choose stopping criterion
epsi = 1e-8;

% parameters shared by the methods
rtol = 1e-4;      % tolerance for the stopping criterion
maxit = 30000;    % maximum number of iterations

%general parameters
verb = 0;  % verbosity

fileid = fopen('bqp-output-TESLA.txt','a');

% BQP problem
bqp_dirs = 'BQP-15000/'; %'BQP-20000/'; 'BQP-25000/';
bqp_dir = "BQP/";
bqp_range = 1:36;


for i = 1:size(bqp_dirs,1)
    for j = bqp_range
        file = strcat('problems/',bqp_dir,bqp_dirs(i,:),'BQP',num2str(j),'.mat');
        load(file)
        fprintf(strcat('BQP',num2str(j)),'\n');
        fprintf('\n');
        fprintf(fileid,'\n');
        fprintf(fileid,strcat('BQP',num2str(j)),'\n');
        fprintf(fileid,'\n');
        
        normsol = norm(x_sol);
        
        mprgpctx.lb = l;
        mprgpctx.ub = u;
        matvec = @(x) MatVetProduct(d,P,x);
        x0 = min(max(x0vect,l),u);
        atol = rtol*norm(matvec(x0)-c);
        
        %% MPRGP
        tic;
        [x1,flg1,k1] = mprgp(matvec,c,x0,atol,maxit,mprgpctx);
        elapsed = toc;
        if verb
            fprintf('Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
            fprintf('MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
                flg1,k1(1),k1(2),k1(3),k1(4),elapsed);
        end
        %
        err_MPRGP = norm(x1 - x_sol)/normsol;
        if verb
            fprintf('Relative error on the solution =%e\n',err_MPRGP);
        end
        fprintf(fileid,'Alg;\t Solved;\t nHess;\t nCG;\t nExp;\t nProp;\t Time\n');
        fprintf(fileid,'MPRGP;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg1,k1(1),k1(2),k1(3),k1(4),elapsed);
        fprintf(fileid,'Relative error on the solution =%e\n',err_MPRGP);
        tic;
        
        SolvedMPRGP = flg1; nHess_MPRGP = k1(1); nCG_MPRGP = k1(2);
        nExp_MPRGP = k1(3); nProp_MPRGP = k1(4); times_MPRGP = elapsed;
        
        if saveout
            save(['output-data-1/outMPRGP_BQP',num2str(j),'.mat'],...
                'SolvedMPRGP', 'nHess_MPRGP', 'nCG_MPRGP','nExp_MPRGP',...
                'nProp_MPRGP','err_MPRGP','times_MPRGP');
        end
        
        %% MPRGPp
        [x2,flg2,k2] = mprgp2(matvec,c,x0,atol,maxit,mprgpctx);
        elapsed = toc;
        if verb
            fprintf('MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
                flg2,k2(1),k2(2),k2(3),k2(4),elapsed);
        end
        err_MPRGPp = norm(x2 - x_sol)/normsol;
        if verb
            fprintf('Relative error on the solution =%e\n',err_MPRGPp);
        end
        fprintf(fileid,'MPRGPp;\t %d;\t %d;\t %d;\t %d;\t %d;\t %g\n',...
            flg2,k2(1),k2(2),k2(3),k2(4),elapsed);
        fprintf(fileid,'Relative error on the solution =%e\n',err_MPRGPp);
        
        SolvedMPRGPp = flg2; nHess_MPRGPp = k2(1); nCG_MPRGPp = k2(2);
        nExp_MPRGPp = k2(3); nProp_MPRGPp = k2(4); times_MPRGPp = elapsed;
        if saveout
            save(['output-data-2/outMPRGPp_BQP',num2str(j),'.mat'],...
                'SolvedMPRGPp', 'nHess_MPRGPp', 'nCG_MPRGPp','nExp_MPRGPp',...
                'nProp_MPRGPp','err_MPRGPp','times_MPRGPp');
        end
        
        %% GP-BoxVABBmin
        [x3,info,f3,step3,err3,vet_ls3,iterati3,times3,nproj3,nhprod3] = ...
            boxABB_min_quad(matvec,c,l,u,x0,rtol,epsi,...
            'MAXIT',maxit,...
            'VERBOSE',0,...
            'M', M_Armijo,...
            'SOL',x_sol,...
            'OPTION',optioniterati,...
            'SAVE',1,...
            'STOPCRITERION',stopcrit,...
            'PARAMETER',[tau,factor,Mrho]);
        iter3 = length(step3);
        if verb
            fprintf('Alg:\t\t\t info\t   nHess\t  nProj\t  nBack  nIter\t  Time\n');
            fprintf('BoxVABBmin - tau=%g %d %d\t   %d\t   %d\t  %d\t %g\n',tau,...
                info,nhprod3,nproj3,length(vet_ls3),iter3,times3(end));
            fprintf('Relative error on the solution =%e\n',err3(end));
        end
        fprintf(fileid,'Alg:\t info\t   nHess\t  nProj\t  nBack  nIter\t  Time\n');
        fprintf(fileid,'BoxVABBmin\t %d  %d\t   %d\t   %d\t  %d\t %g\n',...
            info,nhprod3,nproj3,length(vet_ls3),iter3,times3(end));
        fprintf(fileid,'Relative error on the solution =%e\n',err3(end));
        
        xBoxVABBmin = x3; risBoxVABBmin = f3; stepBoxVABBmin = step3;
        errBoxVABBmin = err3; vet_lsBoxVABBmin = vet_ls3;
        timesBoxVABBmin = times3(end);  nprojBoxVABBmin = nproj3;
        nhprodBoxVABBmin = nhprod3;
        
        if saveout
            save(['output-data-3/outBoxVABBmin_BQP',num2str(j),'.mat'],...
                'tau','xBoxVABBmin', 'risBoxVABBmin','stepBoxVABBmin',...
                'errBoxVABBmin','vet_lsBoxVABBmin','timesBoxVABBmin',...
                'nprojBoxVABBmin','nhprodBoxVABBmin')
        end
        
    end
end
fclose(fileid);
rmpath('mprgp')
rmpath('utils')
rmpath('boxABBmin_solvers')