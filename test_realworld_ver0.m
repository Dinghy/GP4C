%% Run GP4C,GP4CW,LocalEM on the real world data sets.
clear; close all;clc;
addpath(genpath('./'));
% dbstop if warning
% dbstop if error
warning('off','all');
%% Settings
% data generation settings
options.dataReal = 1;           % 1 real data, 0 synthetic data
% training settings
options.maxIter= 1000;          % maximum iteration number 
options.debug = 0;              % 1 debug mode ; 0 do not debug
options.hyper = 1;              % 1 tune hyper ; 0 do not tune
options.hyperOpt = 0;
options.sampleRepeat = 40;      % repeat times for sampling different test data sets.
options.pseudo = 18;            % pseudo inputs in GP
% figure setting
options.figure = 0;             % 1 show figure when generating data; 0 do not show
options.detailPlot = 0;
options.evaApprox = 0;          % 1 see how approximation goes
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood
%% Training and testing
strFileSave = {'RealNausea.mat','RealBladder.mat','RealSkin.mat'};
for idataset = 1:length(strFileSave)
    options.dataType = idataset;           % 1 Nausea, 2 Bladder, 3 Skin
    switch idataset
        case {1,2}
            nPart = 2;
        case 3
            nPart = 4;
    end
    % Initialization
    resFull = cell(nPart,1);   % Different part
    resFull_rev = cell(nPart,1);
    for ip = 1:nPart
        resFull{ip} = cell(2,options.sampleRepeat);
        resFull_rev{ip} = cell(2,options.sampleRepeat);
    end
    % training and testing on each part
    for is = 1:options.sampleRepeat    
        % generate the train-test splitting
        [dataFull,strfileName] = generateRealTest(options);
        
        % obtain the reversed data splitting
        dataFull_rev = reverseSplitting(dataFull);
        
        % process the data and the reversed data
        resFull = runRealworld(dataFull,resFull,nPart,options,is,idataset);
        resFull_rev = runRealworld(dataFull_rev,resFull_rev,nPart,options,is,idataset);
        
        fprintf('\n');
    end
    save(strFileSave{idataset},'resFull','resFull_rev');
end
