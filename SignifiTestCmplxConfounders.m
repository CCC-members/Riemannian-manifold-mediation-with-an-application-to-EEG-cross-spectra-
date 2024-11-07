clear all; close all; clc

addpath(genpath('/home/Carlos/matlab/tensor_toolbox')) %<-- Add the tensor toolbox to the MATLAB path
% cd met; addpath('/home/Carlos/matlab/tensor_toolbox') %<-- [OPTIONAL] Also add the met directory
addpath(genpath('/home/Carlos/matlab/TensorReg')) %<-- Add the toolbox to the Matlab path
% addpath(genpath('/home/Carlos/matlab/TensorReg/private'))
addpath(genpath('/home/Carlos/matlab/mex files'))
addpath(genpath('D:\User1 OneDrive\OneDrive - CCLAB\Thesis\Software\mex files'));
startTime = tic;
sim=1; %Use simulation data (cross and asterisk) or real data
do_boots=1; %Whether to do bootstrap or load the bootstraped data
nboot = 1000;
w=1; % number of confounders

if sim
    n=100;
    % [x,M,y,p,q]=DataGenerator_complex(n); % M out as (p,q,n)
    [x,M,y,p,q]=DataGenerator_complexCovariatesOptimized(n,w+1);
    if isreal(M)
        M_real_im = M;
    else
        M_real_im = [real(M); imag(M)];
    end
    p1 = size(M_real_im,1);
    p=size(M,1);
    q=size(M,2);
    p0 =  2;    % number of confounder
    % confounders = randn(n,p0);
    % x = [x,confounders];
else
    % load('Data78AgeGenderConfounders.mat');
    load('Data78AgeGenderSesConfounders.mat');
    if isreal(M)
        M_real_im = M;
    else
        M_real_im = [real(M); imag(M)];
    end
    p = size(M,1);
    p1 = size(M_real_im,1);
    q = size(M,2);
    n =size(M,3);
    % x=[x,age,gender];
    % x=[x];
    % x=[x,gender,age,ses];
end
numbLambd=30;
Z_design=ones(n,1);
% [A, B, B_all, AIC1, lambdas, A_all, AIC_A, lamdas_A]=Regressions_CovManual(x, M, y, p, q, n, numbLambd); 
[A, B,C_matrices]=Regressions_CovConf(x, M, y, p, q, n, numbLambd);
% A=complex(A(1:p,:),A(p+1:end,:));
B=complex(B(1:p,:),B(p+1:end,:));
AB=A.*B;
if do_boots
    [~,bootsam] = bootstrp(nboot,[],x);
    Aboot=zeros(p,q,nboot);
    Bboot=zeros(p*2,q,nboot);
    parfor i=1:nboot
        x_boot = x(bootsam(:,i),:);
        y_boot = y(bootsam(:,i));
        M_boot = M(:,:,bootsam(:,i));
        [Aboot(:,:,i), Bboot(:,:,i)]=Regressions_CovConf(x_boot, M_boot, y_boot, p, q, n, numbLambd);
    disp(i)
    end
else
    if sim
        load('boot_5000.mat');
    else
        load('RD_boot200.mat');
    end
end
elapsedTime = toc(startTime);
minutes = elapsedTime / 60;
disp(['Elapsed time: ' num2str(minutes) ' minutes.']);

% Aboot=complex(Aboot(1:p,:,:),Aboot(p+1:end,:,:));
Bboot=complex(Bboot(1:p,:,:),Bboot(p+1:end,:,:));
ABboot=Aboot.*Bboot;

num_boots=size(Aboot,3);
AbootStd = zeros(p,q,num_boots);
BbootStd = zeros(p,q,num_boots);
ABbootStd = zeros(p,q,num_boots);

% Find Probabilities of each element 
complex_modA = zeros(p,q);
PA = zeros(p,q);
for i=1:p
    for j=1:q
        complex_modA(i,j) = abs(A(i,j))/std(Aboot(i,j,:));
        PA(i,j)=raylcdf(complex_modA(i,j),1,'upper');
    end
end
complex_modB = zeros(p,q);
PB = zeros(p,q);
for i=1:p
    for j=1:q
        complex_modB(i,j) = abs(B(i,j))/std(Bboot(i,j,:));
        PB(i,j)=raylcdf(complex_modB(i,j),1,'upper');
    end
end
complex_modAB = zeros(p,q);
PAB = zeros(p,q);
for i=1:p
    for j=1:q
        complex_modAB(i,j) = abs(AB(i,j))/std(ABboot(i,j,:));
        PAB(i,j)=raylcdf(complex_modAB(i,j),1,'upper');
    end
end
% autoArrangeFigures()


% Find the thresholds for A 
Q = 0.05;
[pID_A,pN_A] = FDR(reshape(PA,1,[]),Q);
pID_A_abs=abs(pID_A)
% Find the thresholds for B 
[pID_B,pN_B] = FDR(reshape(PB,1,[]),Q);
pID_B_abs=abs(pID_B)
% Find the thresholds for AB 
[pID_AB,pN_AB] = FDR(reshape(PAB,1,[]),Q);
pID_AB_abs=abs(pID_AB)

figure
subplot(2,3,1)
imagesc(abs(A)),colorbar;title('Estimated A');xlabel('Channels');ylabel('Frequency (Hz)'),colormap(hot)
subplot(2,3,2)
imagesc(abs(B)),colorbar;title('Estimated B');xlabel('Channels');ylabel('Frequency (Hz)'),colormap(hot)
subplot(2,3,3)
imagesc(abs(AB)),colorbar;title('Estimated A.*B');xlabel('Channels');ylabel('Frequency (Hz)'),colormap(hot)
subplot(2,3,[4 5 6]); % Create dummy subplots to create space
subplot(2,3,5)
imagesc(abs(C_matrices{1,1})),colorbar;title('Estimated W1');xlabel('a.u.');ylabel('a.u.'),colormap(hot)

figure, 
uPA = PA;
uPA(PA>pID_A) = 0;
uPA(PA<=pID_A) = 1;
subplot(2,2,1)
imagesc(uPA),colorbar;title('Significance of Estimated A');

uPB = PB;
uPB(PB>pID_B) = 0;
uPB(PB<=pID_B) = 1;
subplot(2,2,2)
imagesc(uPB),colorbar;title('Significance of Estimated B');

uPAB = PAB;
uPAB(PAB>pID_AB) = 0;
uPAB(PAB<=pID_AB) = 1;
subplot(2,2,3)
imagesc(uPAB),colorbar;title('Significance of Estimated AB');
