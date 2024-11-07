function [A,B1,C_matrices]=Regressions_CovConf(x,M,y,p,q,n,numbLambd)
%   Regressions_Cov Performs two matrix regressions for a Mediation Model:
%   M = X*A + E, and 
%   y = x*c + M*B + E
%
%   [A,B]=Regressions_Cov(x,M,y,p,q,n,numbLambd)
%   INPUTS:
%   x           stimulus data (vector of 0s and 1s)
%   M           Mediation data (matrix of complex values)
%   y           response data (vector of double values)
%   p           number of rows of M
%   q           number of columns of M
%   numbLambd   number of lambda points to be evaluated to find the model
%               with the minimum AIC
%   OUTPUTS:
%   A           matrix regression coeficient 
%   B           matrix regression coeficient 

%% M = A*x regression
lambdavalue=linspace(-5,5,numbLambd);
lambda=exp(lambdavalue);
DOF=zeros(size(lambda,2),1);
AIC=zeros(size(lambda,2),1);
BIC=zeros(size(lambda,2),1);
lossvalue=zeros(size(lambda,2),1);
% A=zeros(p,q,size(lambda,2),1);
s=size(x,2);

% Center the response and standardize the covariates
Ms = permute(M,[1 3 2]);
Ms = reshape(Ms,[],size(M,2),1);   %convert from 3D to 2D
meanM=mean(M,3);
meanX=mean(x,1);   
stdtempX=std(x,1);
for i=1:n
    Mc((p*i-p+1):(p*i),:)=Ms((p*i-p+1):(p*i),:)-meanM;
    xstd(i,:)=(x(i,:)-meanX)./stdtempX;
end

for i=1:numbLambd
%%% the final estimate of A corresponding to the selected x 
    % [estA(:,:,i),DOF(i)]=RegMatrixCVestimationDOFcov(xstd,Ms,p,q,s,...
    %                                            n,lambda(i));
    [estA(:,:,i),DOF(i)]=RegMatrixCVestimationDOFcovOpt(xstd,Ms,p,q,s,...
                                               n,lambda(i));
    ndf = n*p*q; 
    lossvalue(i) = ErrorLossmatrixdividen(Ms,xstd,estA(:,:,i),s,ndf);
    AIC(i) = lossvalue(i) + 2*DOF(i);
    BIC(i) = lossvalue(i) + log(ndf)*DOF(i);
end
[~,minAIC]=min(AIC);
[~,minBIC]=min(BIC);
A1 = estA(:,:,minAIC);
if(isempty(find(A1)) && minAIC>1)
    A1 = estA(:,:,minAIC-1);
end

for j=1:s
    if s==1
        C_matrices = 0;
    end
    if j==1
        A=A1(:,(1+q*(j-1)):(j*q));
    else
        C_matrices{j-1}=A1(:,(1+q*(j-1)):(j*q));
    end
end

Mc = permute(reshape(Mc.', q,p,[]),[2 1 3]); %convert from 2D to 3D
%% y = M*B regression
numbLambd = 30;
lambdas = zeros(1,numbLambd);
gs = 2/(1+sqrt(3));
B = cell(1,numbLambd);
AIC1 = zeros(1,numbLambd);
BIC1 = zeros(1,numbLambd);
AICn = zeros(1,numbLambd);
BICn = zeros(1,numbLambd);
AII = zeros(1,numbLambd);
BII = zeros(1,numbLambd);
DOF1 = zeros(1,numbLambd);
yhat = zeros(n,1,numbLambd);
RMSE1 = zeros(1,numbLambd);
if isreal(M)
    M_real_im = M;
else
    M_real_im = [real(M);imag(M)];
end
Mt = tensor(M_real_im);
[~,~,stats] = matrix_sparsereg(xstd,Mt,y,inf,'normal','penalty','enet','penparam',1.5);% Eleastic-net penalty (values between 1 and 2)
maxlambda = stats.maxlambda*.95;
for i=1:numbLambd
    if (i==1)
        B0 = [];
    else
        B0 = B{i-1};    % warm start
    end
    lambda = maxlambda*gs^(i-1);
    lambdas(i) = lambda;
    [~,B{i},stats] = matrix_sparsereg(xstd,Mt,y,lambda,'normal','B0',B0,'penalty','enet','penparam',1.5);% Eleastic-net penalty (values between 1 and 2)
    % Collect statistics
    yhat(:,:,i) = stats.yhat;
    AIC1(i) = stats.AIC;
    BIC1(i) = stats.BIC;
    yest = double(ttt(tensor(double(B{i})), Mt, 1:2));
    RMSE1(i) = sqrt(sum((y-yest).^2)/n);
%     AICn(i) = n*log(sum(((y-yest).^2)/n)) + 2*stats.dof;
%     BICn(i) = n*log(sum(((y-yest).^2)/n)) + log(n)*stats.dof;
%     AII(i) = sum((y-yest).^2) + 2*stats.dof;
%     BII(i) = sum((y-yest).^2) + log(n)*stats.dof;
%     DOF1(i) = stats.dof;
end

[~,minAIC1]=min(AIC1);
[~,minBIC1]=min(BIC1);
% loglog(lambdas,AIC1,'-+', lambdas,BIC1, '-o') 
% hold all 
% loglog(lambdas(minAIC1),AIC1(minAIC1),'*b')
% loglog(lambdas(minBIC1),BIC1(minBIC),'*r')
% legend('AIC', 'BIC', 'Location', 'northwest');
% xlabel(' \lambda')
% title('AIC and BIC vs lambda')
B1 = double(B{minAIC1});
end