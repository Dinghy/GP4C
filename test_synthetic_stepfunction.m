%% A simple file to test everything works or not
clear; close all;clc;
addpath(genpath('./'));
% dbstop if warning
% dbstop if error
warning('off','all');
%% settings
% for data generation
options.dataReal = 0;           % 1 real data; 0 synthetic data
options.dataPlot = 0;           % 1 plot when generating data
options.censorAll = 1;          % 1 censor the entire sequence
options.censorNum = 10;         % n number of censored periods when censoring all
options.sampleRepeat = 40;      % n repeat times 
options.censorExperiment = 100; % n how many subjects in 1 repeat time 
options.dataType = 1;           % 1 step function
% training setting
options.hyperOpt = 0;           % how many hyperparameters in a GP we want to optimize
options.debug = 0;              % 1 debug mode; 0 do not debug
options.figure = 0;             % 1 plot the change of estimation during inference
options.detailPlot = 0;         % 1 plot the details of parameters after inference
options.bench = 1;              % 1 run benchmark method
options.benchPos = 15;          % n bench method check points
options.pseudo = 30;            % n number of pseudo inputs GP4C,GP4CW
options.trainMethod = [1,2,3,4,6];  % Test GP3, GP4C(1,0,0.3) and LocalEM
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood

%% Check two different settings of data 
strSave = 'ResultSynA.mat';
options.dataWeight = 0;         % 1 add variations; 0 no variations.
% space to save the result
cRes.True = cell(1,options.sampleRepeat);
dataFull = cell(1,options.sampleRepeat);
for is = 1:options.sampleRepeat
    disp(is);
    % generate the data set first
    dataFull{is} = generateIntCenBook(options);
    cRes.True{is} = testmodel(dataFull{is},[],options);
    fprintf('True MISE %.4f\t %.4f\n',cRes.True{is}.nMISE,cRes.True{is}.nTestLike);
end
%% perform inference on the data set
[sTmp,cModel,vBest] = runSynthetic(options,dataFull);
cRes.Test = sTmp;
%% save the result mat
save(strSave,'cRes','cModel','vBest');

