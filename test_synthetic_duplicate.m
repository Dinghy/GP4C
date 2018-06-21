%% Test file to test the performance of different traing propotions
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
options.censorExperiment = 140; % number of subjects
options.sampleRepeat = 40;       % n repeat times  
options.dataType = 1;           % 1 step function
% training setting
options.hyperOpt = 0;           % set this to 0 
options.debug = 0;              % 1 debug mode; 0 do not debug
options.bench = 1;              % 1 run benchmark method
options.benchPos = 15;          % number of points to be tested in benchmark method
options.pseudo = 30;            % number of pseudo inputs
options.trainMethod = [4,6];    % only test GP4C and LocalEM
% test settings
options.testSample = 3001;      % Sampling point on the whole space for mise and test likelihood
options.testSampTimes = 50;     % Sampling times for MISE and test likelihood
%% Specific testing 
arrDuplicate = [1,0];           % Ratio of training data sets

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
rng(0); 
cRes.Test = cell(1,length(arrDuplicate));
cRes.Duplicate = zeros(length(arrDuplicate),options.sampleRepeat);
for i = 1:length(arrDuplicate)
    np = arrDuplicate(i);
    % change the number of duplicates points
    dataTmp = dataFull;
    for is = 1:options.sampleRepeat
        arrdup = [];
        for j = 1:dataTmp{is}.nTrain
            % round each number with some possibility
            for k = 1:dataTmp{is}.censor{j}.nCInt-1 
                nfloor = floor(dataTmp{is}.censor{j}.vCInt(k,2));
                if rand()<=np && nfloor > dataTmp{is}.censor{j}.vCInt(k,1) && nfloor < dataTmp{is}.censor{j}.vCInt(k+1,2)
                    dataTmp{is}.censor{j}.vCInt(k,2) = nfloor;
                    dataTmp{is}.censor{j}.vCInt(k+1,1) = nfloor;
                end
            end
            % re-form the data
            for k = 1:dataTmp{is}.censor{j}.nCInt
                nBeg = dataTmp{is}.censor{j}.vCInt(k,1);
                nEnd = dataTmp{is}.censor{j}.vCInt(k,2);
                dataTmp{is}.censor{j}.vLen(k) = nEnd-nBeg;
                dataTmp{is}.censor{j}.vCen(k) = sum(dataTmp{is}.train{j}.X >= nBeg & dataTmp{is}.train{j}.X < nEnd); 
            end
            arrdup = [arrdup;dataTmp{is}.censor{j}.vCInt];
        end
        cRes.Duplicate(i,is) = length(unique(arrdup(:,1)));
    end
    % do the training and testing
    sTmp = runSynthetic(options,dataTmp);
    cRes.Test{i} = sTmp;
end
save('SyntheticDuplicate.mat','cRes');