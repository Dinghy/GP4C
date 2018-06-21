%% Plot Figure 4 in the paper
% Illustrate why b=1 is not a good choice.
clear;close all;clc;
addpath(genpath('.\'));

%% Precomputation
load('ResultSynA.mat');
cSave = cell(1,3);
% b = 1
model = cModel{2};
options.testSample = 3001;
vT = linspace(model.Tstart,model.Tend,options.testSample)';
Kmm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.Xm);
Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
KnmInvKmm = Knm/Kmm;
vMean = KnmInvKmm*model.var.Mu;
vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);

cSave{1}.vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
cSave{1}.vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
cSave{1}.vMean = vMean.^2+vVar;

% b = 0.3
model = cModel{4};
Kmm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.Xm);
Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
KnmInvKmm = Knm/Kmm;
vMean = KnmInvKmm*model.var.Mu;
vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);

cSave{2}.vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
cSave{2}.vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
cSave{2}.vMean = vMean.^2+vVar;

% true
fBasic = @(x)7*(x>=0 & x<10)+2*(x>=10 & x<20)+7*(x>=20 & x<30)...
            +2*(x>=30 & x<40)+7*(x>=40 & x<50)+2*(x>=50 & x<60);
cSave{3}.vMean = fBasic(vT);
save('plot\MiddleMat\Figure4.mat','cSave');

%% Plot
close all;
load('plot\MiddleMat\Figure4.mat');load('ResultSynA.mat');model = cModel{2};
fig = figure;set(fig,'Position',[0,0,600,300]);
subplot('Position',[0.09,0.18,0.87,0.70]);hold on;grid on;box on;
% Plot b = 1
ilinewidth = 1.5;
options.testSample = 3001;
vT = linspace(model.Tstart,model.Tend,options.testSample)';

h1 = plot(vT,cSave{1}.vMean,'r','LineWidth',ilinewidth);
plot(vT,cSave{1}.vHigh,'r--','LineWidth',ilinewidth);
plot(vT,cSave{1}.vLow,'r--','LineWidth',ilinewidth);
% plot b = 0.3
h2 = plot(vT,cSave{2}.vMean,'b','LineWidth',ilinewidth);
plot(vT,cSave{2}.vHigh,'b--','LineWidth',ilinewidth);
plot(vT,cSave{2}.vLow,'b--','LineWidth',ilinewidth);
% plot True
h3 = plot(vT,cSave{3}.vMean,'k','LineWidth',ilinewidth);
leg = legend([h1,h2,h3],'GP4C (b=1)','GP4C (b=0.3)','True');
set(leg,'Position',[0.09,0.92,0.87,0.05],'Orientation','Horizontal');
xlabel('Time x');ylabel('Intensity');ylim([0 12]);
set(gca,'FontSize',13);legend boxoff;
printFig(fig,'plotResult\BadB.pdf');
