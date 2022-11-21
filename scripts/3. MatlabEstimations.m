%Spatial estimation, Spatial Economics Analysis Journal R&R

% read data
clear
clc
W1 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\Wq1.txt');
W2 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\Wq2.txt');
W3 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\Wq3.txt');
W4 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\WDk.txt');
W5 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\W1k.txt');
W6 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\W2k.txt');
W7 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\W3k.txt');
W8 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\W4k.txt');
W9 = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\W5k.txt');

A  = readtable('C:\Users\Leonardo\Dropbox\Spatial_Analisis_IC\LM\Stata_CI_Paper3\do files Paper3 CI\EstimacionesLM\Data.xls');

%=====Matrices=====%
W1 = W1{:,2:end};
W2 = W2{:,2:end};
W3 = W3{:,2:end};
W4 = W4{:,2:end};
W5 = W5{:,2:end};
W6 = W6{:,2:end};
W7 = W7{:,2:end};
W8 = W8{:,2:end};
W9 = W9{:,2:end};
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
display('Results Binary Contiguity matrix')
%
% All models include spatial and time period fixed effects (sfe and tfe), (info.)model=3
% model=0 no sfe and tfe, model=1 only sfe, model=2 only tfe
% ----------------------------------------------------------------------------------------
% OLS model without spatial interaction effects
model=3; % sfe + tfe
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('deltaIC','L1_er_','L1_asp_','conflict','dpark','dresguardo','infrastructurepcc', ...
             'humancappcc','induscomerpcc','gaspcc','nontaxincomepcc','royaltiespcc','pserv'); % should be changed if x is changed
prt_reg(results,vnames)

% Additional results
et=ones(T,1);
en=ones(N,1);
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
loglik=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid

% Direct and indirect effects already determined by OLS
% ----------------------------------------------------------------------------------------
% Spatial lag model SAR
info.lflag=0; % required for exact results
info.model=3; % sfe + tfe
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
results=sar_panel_FE(y,x,W,T,info); 
% Print out coefficient estimates
vnames=strvcat('deltaIC','L1_er_','L1_asp_','conflict','dpark','dresguardo','infrastructurepcc',...
              'humancappcc','induscomerpcc','gaspcc','nontaxincomepcc','royaltiespcc','pserv'); % variable names               
prt_sp(results,vnames,1);
% Print out direct and indirect effects
spat_model=0; % 0 for SAR, 1 for SDM/GNS model
direct_indirect_effects_estimates(results,W,spat_model);


% ----------------------------------------------------------------------------------------
% Spatial error model SEM
info.lflag=0; % required for exact results
info.model=3; % sfe + tfe
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
results=sem_panel_FE(y,x,W,T,info); 
% Print out coefficient estimates
vnames=strvcat('deltaIC','L1_er_','L1_asp_','conflict','dpark','dresguardo','infrastructurepcc',...
    'humancappcc','induscomerpcc','gaspcc','nontaxincomepcc','royaltiespcc','pserv'); % variable names
prt_sp(results,vnames,1);
% Direct and indirect effects already determined by SEM
% ----------------------------------------------------------------------------------------
% SLX Model  ... NA
% Direct and indirect effects already determined by SLX
% ----------------------------------------------------------------------------------------
% Spatial Durbin model SDM
info.lflag=0; % required for exact results
info.model=3; % sfe + tfe
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
info.bc=1; % bias correction of Leen and Yu on
results=sar_panel_FE(y,[x wx],W,T,info); 
% Print out coefficient estimates
vnames=strvcat('deltaIC','L1_er_','L1_asp_','conflict','dpark','dresguardo','infrastructurepcc',...
             'humancappcc','induscomerpcc','gaspcc','nontaxincomepcc','royaltiespcc','pserv',...
              'W*L1_er_','W*L1_asp_','W*conflict','W*dpark','W*dresguardo','W*infrastructurepcc',...
              'W*humancappcc','W*induscomerpcc','W*gaspcc','W*nontaxincomepcc','W*royaltiespcc','W*pserv'); % variable names
prt_sp(results,vnames,1);
% Print out direct and indirect effects
spat_model=1;
direct_indirect_effects_estimates(results,W,spat_model);
% ----------------------------------------------------------------------------------------
% Spatial Durbin error model
info.lflag=0; % required for exact results
info.model=3; % sfe + tfe
info.fe=0; % Do not print intercept and fixed effects; use info.fe=1 to turn on
info.bc=1; % bias correction of Leen and Yu on
results=sem_panel_FE(y,[x wx],W,T,info); 
% Print out coefficient estimates
vnames=strvcat('deltaIC','L1_er_','L1_asp_','conflict','dpark','dresguardo','infrastructurepcc',...
             'humancappcc','induscomerpcc','gaspcc','nontaxincomepcc','royaltiespcc','pserv',...
              'W*L1_er_','W*L1_asp_','W*conflict','W*dpark','W*dresguardo','W*infrastructurepcc',...
              'W*humancappcc','W*induscomerpcc','W*gaspcc','W*nontaxincomepcc','W*royaltiespcc','W*pserv'); % variable names
prt_sp(results,vnames,1);
% Direct and indirect effects already determined by SDEM
% ----------------------------------------------------------------------------------------
% General nesting spatial model GNS
% ----------------------------------------------------------------------------------------

% Bayesian approach to determine best combination of model and W matrx
% beta prior
display('Bayaesian results per W matrix')
prior.c = 1.01;
prior.d = 1.01;
opt.tol = 1e-3; opt.disp = 0;
incr = 0.001;
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
ted = 1; % set ted=0 for model with spatial fixed effects without time dummies, % set ted=1 for model with spatial and time period fixed effects
model=3;
[yf,xf,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

W=W1; % Binary contiguity matrix
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[result.lmarginal];

W=W2;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

W=W3;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]


W=W4;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]

W=W5;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]


W=W6;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]



W=W7;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]



W=W8;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]

W=W9;
[n,junk] = size(W);
lambda = eigs(sparse(W),speye(n),1,'SR',opt); 
rmin = real(1/lambda) + incr;   
rmax = 1.0 - incr;
[rmin rmax]
out=lndetfull(sparse(W),rmin,rmax);
result = log_marginal_panelprob(yf,xf,W,N,T,model,out,prior,rmin,rmax,incr); 
in.cnames = strvcat('W-miles log marginals','model probs');
in.rnames = strvcat('Models','sar','sdm','sem','sdem');
mprint([result.lmarginal result.probs],in);
lmarginal=[lmarginal;result.lmarginal];

nmodels = length(lmarginal);
adj = max(lmarginal(:,1));
madj = matsub(lmarginal,adj);
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs = [matdiv(xx,psum)];
display('Overall results models and W; 4 rows (sar,sdm,sem,sdem) for each W')
[lmarginal/1000 probs]

