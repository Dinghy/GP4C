%% Plot Figure 11 in the appendix
% Comparison of the performance of GP3, GP4C, and LocalEM when varying the
% number of duplicates

clear;close all;clc;
addpath(genpath('./.'));
load('SyntheticDuplicate.mat');
vX = [1,0.8,0.6,0.4,0.2,0];
strleg = {'GP3','GP4C','LocalEM'};
%%
fig = figure;
set(fig,'Position',[100,100,600,350]);

options.errbar = 1;
plotEbar(vX,cRes.Duplicate',fig,[0,0,1],[0.12,0.15,0.38,0.78],options);
hold on;grid on;box on;set(gca,'FontSize',13);
axis([-0.1,1.1,0,700]);
set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1]);
title('Duplicate Number');
xlabel('Probability $$p_0$$','Interpreter','latex');
ylabel('The number $$\bar{N}$$','Interpreter','latex');

vmethod = [4,6];
mcolor = [1,0,0;0,0,1];
for i = 1:length(vmethod)
    mTime = zeros(length(vX),30);
    for j = 1:length(vX)
        mTime(j,:) = cellfun(@(x)x.nTime,cRes.Test{j}(vmethod(i),:));
    end
    plotEbar(vX,mTime',fig,mcolor(i,:),[0.6,0.15,0.38,0.78],options);
end
hold on;grid on;box on;set(gca,'FontSize',13);
axis([-0.1,1.1,0,60]);
set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1]);
title('Computation Time');
xlabel('Probability $$p_0$$','Interpreter','latex');ylabel('time (seconds)');
legend('GP4C','LocalEM');
printFig(fig,'result_plot\SynDuplicate'); % Figure 8 in the appendix