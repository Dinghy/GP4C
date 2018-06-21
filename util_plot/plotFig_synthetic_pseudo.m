%% Plot Figure 6 in the paper
% illustrate the performance of different pseudo inputs

clear;close all;clc;
addpath(genpath('./.'));

%% File Ratio
fig = figure;set(fig,'Position',[0,0,1000,350]);
load('SyntheticPseudo.mat');
vX = 15:5:50;
mColor = [1,0.8,0;1,0,0;0,0,1];
vPlotMethod = [6,1,4];
strleg = {'LocalEM','GP3','GP4C'};
%% test likelihood
vPos = [0.07,0.24,0.25,0.66];
arrTrue = cellfun(@(x)x.nTestLike,cRes.True)';
options.errbar = 0;options.marker = 0;
% plotEbar(vX,arrTrue,fig,[1,0,1],vPos,options);
for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    if i == max(vPlotMethod)
        vY = cellfun(@(x)x.nTestLike,cRes.Test{1}(i,:));
        options.errbar = 1;options.marker = 0;
        plotEbar(vX,vY,fig,mColor(iPos,:),vPos,options);
    else
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
end
subplot('Position',vPos);box on;grid on;
title('Test Likelihood');xlabel('Pseudo Inputs Number');
leg = legend(strleg([1,2,3]),'Orientation','horizontal');
set(leg,'Position',[0,0,1,0.05]);
legend boxoff;
xlim([10 55]);ylim([-1410 -1340]);
set(gca,'FontSize',13);

%% MISE
vPos = [0.4,0.24,0.25,0.66];
for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    if i == max(vPlotMethod)
        vY = cellfun(@(x)x.nMISE,cRes.Test{1}(i,:));
        options.errbar = 1;options.marker = 0;
        plotEbar(vX,vY,fig,mColor(iPos,:),vPos,options);
    else
        mY = zeros(size(cRes.Test{1},2),length(vX));
        for j = 1:length(vX)
            mY(:,j) = cellfun(@(x)x.nMISE,cRes.Test{j}(i,:));
        end
        options.errbar = 1;options.marker = 1;
        plotEbar(vX,mY,fig,mColor(iPos,:),vPos,options);
    end
end
subplot('Position',vPos);box on;grid on;
title('MISE');xlabel('Pseudo Inputs Number');
xlim([10 55]);ylim([10 70]);
set(gca,'FontSize',13);

%% Computation Time
vPos = [0.72,0.24,0.25,0.66];

for iPos = 1:length(vPlotMethod)
    i = vPlotMethod(iPos);
    if i == max(vPlotMethod)
        vY = cellfun(@(x)x.nTime,cRes.Test{1}(i,:));
        options.errbar = 1;options.marker = 0;
        plotEbar(vX,vY,fig,mColor(iPos,:),vPos,options);
    else
        mY = zeros(size(cRes.Test{1},2),length(vX));
        for j = 1:length(vX)
            mY(:,j) = cellfun(@(x)x.nTime,cRes.Test{j}(i,:));
        end
        options.errbar = 1;options.marker = 1;
        plotEbar(vX,mY,fig,mColor(iPos,:),vPos,options);
    end
end
subplot('Position',vPos);box on;grid on;
title('Computation Time');xlabel('Pseudo Inputs Number');
set(gca,'FontSize',13);
xlim([10 55]);
printFig(fig,'result_plot\SynAPseudo');