%% Test file to test the performance of different pseudo inputs
clear; close all;clc;
addpath(genpath('./'));
% dbstop if warning
% dbstop if error
warning('off','all');
%% setting 
% for data generation
options.dataReal = 0;           % 1 real data; 0 synthetic data
options.dataPlot = 0;           % 1 plot when generating data
options.censorAll = 1;          % 1 censor the entire sequence
options.censorNum = 10;         % number of censored periods when censoring all
options.dataWeight = 0;         % add extra weight to each subject when generating data
options.censorExperiment = 100;
options.sampleRepeat = 50;      % n repeat times  
options.dataType = 1;           % 1 step function
% training setting
options.detailPlot = 0;
options.figure = 0;
options.hyperOpt = 0;
options.debug = 0;              % 1 debug mode; 0 do not debug
options.bench = 1;              % 1 run benchmark method
options.benchPos = 15;
options.trainMethod = [1,4,6];  % only test GP3, GP4C and LocalEM
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood
%% Specific testing 
arrPseudo = 15:5:50;   % pseudo inputs number

%% generate the data set first
cRes.True = cell(1,options.sampleRepeat);
dataFull = cell(1,options.sampleRepeat);
for is = 1:options.sampleRepeat
    disp(is);
    dataFull{is} = generateIntCenBook(options);
    cRes.True{is} = testmodel(dataFull{is},[],options);
    fprintf('True MISE %.4f\t %.4f\n',cRes.True{is}.nMISE,cRes.True{is}.nTestLike);
end
%% Train and Testing
cRes.Test = cell(1,length(arrPseudo));
for i = 1:length(arrPseudo)
    if i == 1
        options.bench = 1;
    else
        options.bench = 0;
    end
    options.pseudo = arrPseudo(i);
    sTmp = runSynthetic(options,dataFull);
    cRes.Test{i} = sTmp;
end
% save('SyntheticPseudo.mat','cRes');