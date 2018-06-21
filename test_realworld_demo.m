%% Run GP4C,GP4CW,LocalEM on the real world data sets.
clear; close all;clc;
addpath(genpath('./'));
% dbstop if warning
% dbstop if error
warning('off','all');
%% Run three data sets in a row
% data generation settings
options.dataReal = 1;           % 1 real data, 0 synthetic data
% training settings
options.maxIter= 1000;          % maximum iteration number 
options.debug = 0;              % 1 debug mode ; 0 do not debug
options.hyper = 1;              % 1 tune hyper ; 0 do not tune
options.pseudo = 18;            % pseudo inputs in GP
options.hyperOpt = 0;
% figure setting
options.figure = 0;             % 1 show figure when generating data; 0 do not show
options.detailPlot = 0;
options.sampleRepeat = 1;       % repeat times for sampling different test data sets.
options.evaApprox = 0;          % 1 see how approximation goes
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood
%% Training and testing
strFileSave = {'RealNauseaDemo.mat','RealBladderDemo.mat'};
idataset = 2;
options.dataType = 2;           % 1 Nausea, 2 Bladder, 3 Skin
nPart = 2;

% Initialization
resFull = cell(nPart,1);   % Different part
resModel = cell(nPart,3);  % save a model for each model for plotting
for ip = 1:nPart
    resFull{ip} = cell(2,options.sampleRepeat);
end
% training and testing on each part
for is = 1:options.sampleRepeat
    [dataFull,strfileName] = generateRealTest(options);
    for ip = 1:1
        fprintf('Dataset %d\t Part %d\t Sample %d\n',idataset,ip,is);
        data = dataFull{ip};
        data.nU = data.nTrain;
        %% censored data with GP
        nMethod = 1;
        options.inferType = 3;
        % assume the same intensity function
        options.wModel = 0;
        [model,nTime] = varTrainIntCen(data,options);
        if is == 1
            resModel{ip,nMethod} = model;
        end
        res = testmodel(data,model,options);
        res.nTime = nTime;
        displayResult(res);
        resFull{ip}{nMethod,is} = res;
        nMethod = nMethod + 1;
        %% censored data with GP and adding additional weight
        options.wModel = 1;
        options.wLow = 1e-6;
        options.inferType = 3;
        [model,nTime] = varTrainIntCen(data,options);
        if is == 1
            resModel{ip,nMethod} = model;
        end
        res = testmodel(data,model,options);
        res.nTime = nTime;
        displayResult(res);
        resFull{ip}{nMethod,is} = res;
        nMethod = nMethod+1;
        options.wModel = 0;
        %% Benchmark -- localEM
        [model,nTime] = localEM(data,options);
        res = testmodel(data,model,options);
        if is == 1
            resModel{ip,nMethod} = model;
        end
        res.nTime = nTime;
        displayResult(res);
        resFull{ip}{nMethod,is} = res;
    end
    fprintf('\n');
end
save(strFileSave{idataset},'resFull','resModel','data','options');

