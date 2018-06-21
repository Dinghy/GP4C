%% Test synthetic data sets with a intensity function drawn from a Gaussian process
clear; close all;clc;
addpath(genpath('./'));
% dbstop if warning
% dbstop if error
warning('off','all');
%% settings
% for data generation
options.dataReal = 0;           % 1 real data; 0 synthetic data
options.dataPlot = 1;           % 1 plot when generating data
options.censorAll = 1;          % 1 censor the entire sequence
options.censorNum = 10;         % number of censored periods when censoring all
options.censorExperiment = 100;
options.sampleRepeat = 40;      % n repeat times  
options.dataType = 1;           % 1 step function
% training setting
options.hyperOpt = 0;
options.debug = 0;              % 1 debug mode; 0 do not debug
options.figure = 0;
options.detailPlot = 0;
options.bench = 1;              % 1 run benchmark method
options.benchPos = 15;          % bench method check points
options.pseudo = 30;
options.trainMethod = [1,2,3,4,6];        % training method 
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood
% save the result according to the function
%% generate a general function
strList={'A','B'};
for i = 1:2
    vT = linspace(0,60,options.testSample)';
    scov = {'covSum',{'covSEard','covNoise'}};
    Knn = feval(scov{:},[1.0136;-0.2551;-6.9078],vT);
    options.f = chol(Knn,'lower')*randn(size(vT));
    options.f = options.f.^2;
    if options.figure == 1
        figure;plot(vT,options.f);
    end
    save(['sampleGP',strList{i},'.mat'],'options');
    %% Check two different settings of data
    options.dataWeight = 0;         % add variations to the intensity functions in the data sets.
    % generate the data set first
    cRes.True = cell(1,options.sampleRepeat);
    dataFull = cell(1,options.sampleRepeat);
    for is = 1:options.sampleRepeat
        disp(is);
        % load('dataSave1.mat');
        dataFull{is} = generateIntCenBook(options);
        cRes.True{is} = testmodel(dataFull{is},[],options);
        fprintf('True MISE %.4f\t %.4f\n',cRes.True{is}.nMISE,cRes.True{is}.nTestLike);
    end
    % iterate over the number of pseudo inputs
    [sTmp,cModel,vBest] = runSynthetic(options,dataFull);
    cRes.Test = sTmp;
    save(['ResultSynGP',strList{i},'.mat'],'cRes','cModel','vBest');
end