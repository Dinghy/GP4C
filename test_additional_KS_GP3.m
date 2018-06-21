% Comparing Kernel Smoothing method and GP3 method.
clear;close all;clc;
addpath(genpath('./'));
%% options
options.dataType = 2;  % 1 simple line ; 2 mixture curve;
options.dataPlot = 0;  % 1 plot
options.censorNum = 5; 
options.endCorrect = 1;
options.maxIter= 1000;
options.figure = 0; % 1 show figure; 0 do not show
options.debug = 0;  % 1 debug mode ; 0 do not debug
options.hyper = 1;  % 1 tune hyper ; 0 do not tune
options.ratio = 0.1; 
options.sampleRepeat = 30;  % repeat times for sampling
options.pseudo = 18; % pseudo inputs in GP
options.nCurU = 1;
%% test sampling
options.dataPlot = 1;
data = generateIntCen(options);
options.dataPlot = 0;
%% training
res.KSerror = zeros(options.sampleRepeat,1);
res.GPerror = zeros(options.sampleRepeat,1);
for is = 1:options.sampleRepeat
    disp(is);
    data = generateIntCen(options);
    if is == 1
        options.figure = 1; % 1 show figure; 0 do not show
    else
        options.figure = 0; % 1 show figure; 0 do not show
    end
    % Kernel Smoothing
    modelKS = kernelSmoothing(data,options);
    res.KSerror(is) = testMSE(data,modelKS);
    % GPPP
    [modelGP,nTime] = varTrainGP3(data,options);
    res.GPerror(is) = testMSE(data,modelGP);
end
figure;
boxplot([res.KSerror,res.GPerror]);