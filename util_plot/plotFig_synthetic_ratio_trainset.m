%% Plot Figure 8 in the appendix
% Comparison of the performance of GP3, GP4C, and LocalEM when varying the
% ration of training objects.

clear;close all;clc;
addpath(genpath('./.'));
%% File Ratio
fig = figure;set(fig,'Position',[0,0,1000,300]);
load('SyntheticRatio.mat');
vX = 0.2:0.1:1;
mColor = [1,0,0;0,0,1;1,0.8,0];
vPlotMethod = [1,4,6];
strleg = {'GP3','GP4C','LocalEM'};
%% test likelihood
vPos = [0.07,0.24,0.25,0.66];
arrTrue = cellfun(@(x)x.nTestLike,cRes.True)';
options.errbar = 0;options.marker = 0;

for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    mY = zeros(size(cRes.Test{1},2),length(vX)); 
    for j = 1:length(vX)
        if isfield(cRes.Test{j}{i,1},'nSampLike')
            mY(:,j) = cellfun(@(x)x.nSampLike,cRes.Test{j}(i,:));
        else
            mY(:,j) = cellfun(@(x)x.nTestLike,cRes.Test{j}(i,:));
        end
    end
    options.errbar = 1;options.marker = 1;
    plotEbar(vX,mY,fig,mColor(iPos,:),vPos,options);
end
subplot('Position',vPos);box on;grid on;
title('Test Likelihood');xlabel('Training Subjects Ratio');
set(gca,'FontSize',13);
leg = legend(strleg([1,2,3]),'Orientation','horizontal');
set(leg,'Position',[0,0,1,0.05]);
xlim([0.1 1.1]);ylim([-1440 -1340]);
legend boxoff ;

%% MISE
vPos = [0.4,0.24,0.25,0.66];
for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    mY = zeros(size(cRes.Test{1},2),length(vX)); 
    for j = 1:length(vX)
    	mY(:,j) = cellfun(@(x)x.nMISE,cRes.Test{j}(i,:));
    end
    options.errbar = 1;options.marker = 1;
    plotEbar(vX,mY,fig,mColor(iPos,:),vPos,options);
end
subplot('Position',vPos);box on;grid on;
xlim([0.1 1.1]);ylim([20 80]);
title('MISE');xlabel('Training Subjects Ratio');
set(gca,'FontSize',13);

%% Computation Time
vPos = [0.72,0.24,0.25,0.66];

for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    mY = zeros(size(cRes.Test{1},2),length(vX)); 
    for j = 1:length(vX)
    	mY(:,j) = cellfun(@(x)x.nTime,cRes.Test{j}(i,:));
    end
    options.errbar = 1;options.marker = 1;
    plotEbar(vX,mY,fig,mColor(iPos,:),vPos,options);
end
subplot('Position',vPos);box on;grid on;
title('Computation Time');xlabel('Training Subjects Ratio');
xlim([0.1 1.1]);ylim([0 40]);
set(gca,'FontSize',13);

printFig(fig,'result_plot\SynARatio'); % Figure 8 in the appendix
