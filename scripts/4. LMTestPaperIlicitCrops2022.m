%Lagrange Multiplier Tests%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by: J.Paul Elhorst summer 2010
% University of Groningen
% Department of Economics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.
%
% Elhorst JP (2010) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.
%
%Spatial estimation, Spatial Economics Analysis Journal R&R

% read data
clear
clc
W1 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\Wq2.txt');
A  = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\do files Paper3 CI\EstimacionesLM\Data.xls');

%=====Matrices=====%
W1 = W1{:,2:end};
%=====Matrices=====%

%=====Data=====%
A  = A{:,:};
%=====Data=====%

% set dimensions
mi = min(A(1:end,2));
mx = max(A(1:end,2));

T=mx-mi+1; % number of time periods
N=size(A,1)/T; % number of regions

% row-normalize W
%W=normw(W1); % row-normalize binary contiguity matrix <L>: Matrix already normalized
W=W1;
% define dependent and independent variables
y=A(:,[3]); % column number in the data matrix that corresponds to the dependent variable
x=A(:,[4:end]); % column numbers in the data matrix that correspond to the independent variables
% Create WX variables
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=W*x(t1:t2,:);
end

% determine size of the x matrix
[nobs K]=size(x);
xconstant=ones(N*T,1);
display('Results Binary Contiguity matrix')

% ----------------------------------------------------------------------------------------

% ols estimation 
results=ols(y,[xconstant x]);
vnames=strvcat('logcit','intercept','logp','logy');
prt_reg(results,vnames,1);
sige=results.sige*((nobs-K)/nobs);
loglikols=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid

% The (robust)LM tests developed by Elhorst

LMsarsem_panel(results,W,y,[xconstant x]); % (Robust) LM tests

% The lm tests developed by Donald Lacombe
% see http://www.rri.wvu.edu/lacombe/~lacombe.htm

lm1=lmlag_panel(y,[xconstant x],W);
%prt_tests(lm1);

lm2=lmerror_panel(y,[xconstant x],W);
%prt_tests(lm2);

lm3=lmlag_robust_panel(y,[xconstant x],W);
%prt_tests(lm3);

lm4=lmerror_robust_panel(y,[xconstant x],W);
%prt_tests(lm4);

% ----------------------------------------------------------------------------------------
% spatial fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=1;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
sfe=meanny-meannx*results.beta; % including the constant term
yme = y - mean(y);
et=ones(T,1);
error=y-kron(et,sfe)-x*results.beta;
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
logliksfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests
%%
lm1=lmlag_panel(ywith,xwith,W);
%prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
%prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
%prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
%prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% time-period fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=2;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
tfe=meanty-meantx*results.beta; % including the constant term
yme = y - mean(y);
en=ones(N,1);
error=y-kron(tfe,en)-x*results.beta;
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
logliktfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests

lm1=lmlag_panel(ywith,xwith,W);
%prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
%prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
%prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
%prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% spatial and time period fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=3;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
intercept=mean(y)-mean(x)*results.beta; 
sfe=meanny-meannx*results.beta-kron(en,intercept);
tfe=meanty-meantx*results.beta-kron(et,intercept);
yme = y - mean(y);
ent=ones(N*T,1);
error=y-kron(tfe,en)-kron(et,sfe)-x*results.beta-kron(ent,intercept);
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
loglikstfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid

LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests

lm1=lmlag_panel(ywith,xwith,W);
%prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
%prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
%prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
%prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% Tests for the joint significance of spatial and/or time-period fixed effects
LR=-2*(logliktfe-loglikstfe);
dof=N;
probability=1-chis_prb(LR,dof);
% Note: probability > 0.05 implies rejection of spatial fixed effects
fprintf(1,'LR-test joint significance spatial fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);
LR=-2*(logliksfe-loglikstfe);
dof=T;
probability=1-chis_prb(LR,dof);
% Note: probability > 0.05 implies rejection of spatial fixed effects
fprintf(1,'LR-test joint significance time-periode fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);
