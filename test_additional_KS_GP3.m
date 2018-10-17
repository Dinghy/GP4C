% Comparing Kernel Smoothing method and GP3 method.
clear;close all;clc;
addpath(genpath('./'));
%% options
options.dataType = 6;  
options.dataPlot = 0;  % 1 plot
options.censorNum = 5; 
options.endCorrect = 1;
options.maxIter= 1000;
options.figure = 0; % 1 show figure; 0 do not show
options.debug = 0;  % 1 debug mode ; 0 do not debug
options.hyper = 1;  % 1 tune hyper ; 0 do not tune
options.ratio = 0.1; 
options.sampleRepeat = 1;  % repeat times for sampling
options.pseudo = 30; % pseudo inputs in GP
options.nCurU = 1;

% added to pass the new scripts 2018-10-10
options.censorExperiment = 2;
options.dataWeight = 0;
options.wModel = 0;
options.hyperOpt = 1;
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood

%% test sampling
% options.dataPlot = 1;
% data = generateIntCenBook(options);
% options.dataPlot = 0;
%% training
res.KSerror = zeros(options.sampleRepeat,1);
res.GPerror = zeros(options.sampleRepeat,1);
for is = 1:options.sampleRepeat
    disp(is);
    data = generateIntCenBook(options);
    if is == 1
        options.figure = 1; % 1 show figure; 0 do not show
    else
        options.figure = 0; % 1 show figure; 0 do not show
    end
    % Kernel Smoothing
    modelKS = kernelSmoothing(data,options);
    [res.KSerror(is),a] = testMSE(data,modelKS,options);
    % GPPP
    [modelGP,nTime] = varTrainGP3Book(data,options);
    [res.GPerror(is),b] = testMSE(data,modelGP,options);
    disp([res.KSerror(is),res.GPerror(is),a,b]);
    if res.KSerror(is) > res.GPerror(is) || is == 1
        data_save = data;
        modelKS_save = modelKS;
        modelGP_save = modelGP;
    end
end
%%
data = data_save;
model = modelGP_save;
vT = model.plotXm;
Kmm = feval(model.cov{:},model.GP.logtheta,model.Xm);
Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
KnmInvKmm = Knm/Kmm;
vMean = KnmInvKmm*model.var.Mu;
InvKmmSigma = Kmm\model.var.Sigma;
vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+...
    sum(KnmInvKmm.*(Knm*(InvKmmSigma-eye(model.m))),2);
vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
vMean = vMean.^2+vVar;

fig = figure;
ilinewidth = 1.5;
set(fig,'Position',[0,0,950,500]);

subplot('Position',[0.06,0.2,0.4,0.68]);set(gca,'FontSize',11);
grid on;box on;hold on;
l1 = plot(vT,data.fBasic(vT),'k','LineWidth',ilinewidth);
l2 = plot(vT,modelKS_save.flam(vT),'r','LineWidth',ilinewidth);
l3 = plot(vT,vMean,'b','LineWidth',ilinewidth);
plot(vT,vHigh,'b--','LineWidth',ilinewidth);
plot(vT,vLow,'b--','LineWidth',ilinewidth);
xlim([data.XBeg,data.XEnd]);
ylim([0,max(data.fBasic(vT))]);

set(gca,'XTickLabel',[]); ylabel('Intensity');
legend([l1,l2,l3],{'True','Kernel Smoothing','Variational Inference'},'Location','best');
title('Comparison of the Intensity Functions');

% plot points the bottom plot
subplot('Position',[0.06,0.12,0.4,0.05]);set(gca,'FontSize',11);
for ibar=1:length(data.train{1}.X)
    iX=data.train{1}.X(ibar);
    line([iX,iX],[-1,1],'Color','k','LineStyle','-');
end
set(gca,'YTickLabel',[]);
ylabel('Data');
xlim([data.XBeg,data.XEnd]);
xlabel('Time');

subplot('Position',[0.58,0.05,0.4,0.82]);set(gca,'FontSize',11);
boxplot([res.KSerror,res.GPerror],{'Kernel Smoothing','Variational Inference'});

xlabel('Method');ylabel('ISE');
title('Comparison of the ISE');
%%
printFig(fig,'result_plot\KS_GP3'); % save the figure.