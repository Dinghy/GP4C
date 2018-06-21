%% Plot Figure 5 in the paper.
% Inferred intensity function by LocalEM and GP4C method on the synthetic
% datasets.

clear;close all;clc;
addpath(genpath('.\.'));

% load('ResultSynA.mat');
% fprintf('Result for Synthetic A data set\n');
% res = cRes.Test;
% printStat([],res);
%% Load GP-based functions
clc;ilinewidth = 1.5;
fprintf('\n');
strFile = {'ResultSynGPA.mat','ResultSynGPB.mat'};
strTrue = {'sampleGPA.mat','sampleGPB.mat'};
cPre = cell(1,length(strFile));
for i = 1:length(strFile)
    load(strTrue{i});
    load(strFile{i});
    
    fprintf('Result for Synthetic GP data set\n');
    vT = linspace(0,60,options.testSample)';
    cPre{i}.vT = vT;
    cPre{i}.fTrue = options.f;
    model = cModel{6};
    cPre{i}.fLocalEM = model.pre;
    model = cModel{4};
    Kmm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.Xm);
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
    KnmInvKmm = Knm/Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);
    cPre{i}.vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
    cPre{i}.vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
    cPre{i}.vMean = vMean.^2+vVar;
end
save('util_plot\MiddleMat\Syn.mat','cPre','vT','options');
%% plot
close all;
load('util_plot\MiddleMat\Syn.mat','cPre');
fig = figure;
set(fig,'Position',[0,0,600,400]);
strTitle = {'Synthetic B','Synthetic C'};
for i = 1:2
    subplot('Position',[0.09,0.13+(i-1)*0.40,0.85,0.33]);hold on;grid on;box on;
    % true
    h1 = plot(vT,cPre{i}.fTrue,'k','LineWidth',ilinewidth);
    % benchmark
    model = cModel{6};
    h2 = plot(vT,cPre{i}.fLocalEM,'b','LineWidth',ilinewidth);
    % b = 0.3
    h3 = plot(vT,cPre{i}.vMean,'r','LineWidth',ilinewidth);
    plot(vT,cPre{i}.vHigh,'r--','LineWidth',ilinewidth);
    plot(vT,cPre{i}.vLow,'r--','LineWidth',ilinewidth);
    if i== 1
        xlabel('Time x');ylim([0 8]);
        text(1,7,strTitle{i},'FontSize',13);
    else
        text(1,5,strTitle{i},'FontSize',13);
    end
    ylabel('Intensity');set(gca,'FontSize',13);
    
end

leg = legend([h1,h2,h3],'True','LocalEM','GP4C(b=0.3)','Orientation','Horizontal');
set(leg,'Position',[0.09,0.91,0.85,0.05]);legend boxoff;
printFig(fig,'result_plot\SynB.pdf'); % Figure 5