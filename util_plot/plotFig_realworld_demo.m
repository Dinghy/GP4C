%% Plot figures for Figure 1,2,7,9,10 in the paper.
% 1.  Illustrate the panel count data on the patients (Bladder,train)
% 2.  Inferred intensity function by LoalEM and GP4C
% 7.  Illustrate the panel count data on the patients (Bladder,testing)
% 9.  Inferred intensity function by LoalEM and GP4CW (Appendix)
% 10. Illustrate the panel count data on the patients (Appendix,Bladder,testing)
clear;close all;clc;
addpath(genpath('.\.'));
%% GP4C
load('RealBladderDemo.mat');
resA = resFull{1}{1};resA.model = resModel{1,1};
resB = resFull{1}{3};resB.model = resModel{1,3};
options.plotDemo = 1;options.dataWeight = 0;
[fig1,fig2,fig3] = plotDataModel(data,resA,resB,options);
printFig(fig1,'result_plot\DemoPanel');     % Figure 1
printFig(fig2,'result_plot\TestBladder');   % Figure 7
printFig(fig3,'result_plot\DemoIntensity'); % Figure 2
%% GP4CW
resA = resFull{1}{2};resA.model = resModel{1,2};
resB = resFull{1}{3};resB.model = resModel{1,3};
options.plotDemo = 1;options.dataWeight = 1;
[fig1,fig2,fig3] = plotDataModel(data,resA,resB,options);
printFig(fig2,'result_plot\TestBladderW');   % Figure 9
printFig(fig3,'result_plot\DemoIntensityW'); % Figure 10