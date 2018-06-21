function [resTest,cModel,vBest] = runSynthetic(options,dataFull)
% RUN SYNTHETIC DATA EXPERIMENT GIVEN THE DATA SET, 
% Please do not run real world data sets!!!
%
% Input :  datafull -- contains the data sets, if missing, the system will
%                      generate one synthetic data set
%           options -- training and testing options
% Output:   resTest -- cell which contains the result: Size: 4*options.sampleRepeat
%                   field: nTime    -- running time 
%                          kernel   -- exists for benchmark, bandwidth
%                          (others) -- see testmodel.m
%           cModel -- store the best model in terms of MISE during all trials


%% setting
% data generation settings
if nargin == 1  % no data
    options.dataReal = 0;           % 1 real data; 0 synthetic data
    options.dataPlot = 0;           % 1 plot when generating data
    options.censorAll = 1;          % 1 censor the entire sequence
    options.censorNum = 10;         % number of censored periods when censoring all
    % number of experiments
    if options.dataType == 5
        options.censorExperiment = 20;   % number of sequences in one training
    else
        options.censorExperiment = 100;  % number of sequences in one training
    end
    options.sampleRepeat = 1;
    options.dataWeight = 0;              % 1 add extra weight when generating the data   
    dataFull = cell(1,options.sampleRepeat);
    for iSampleTimes = 1:options.sampleRepeat
        dataFull{iSampleTimes} = generateIntCenBook(options);
    end
end
% figure settings
if ~isfield(options,'figure')
    options.figure = 0;
end
% training settings
options.maxIter= 1000;                  % maximum iteration number 
options.hyper = 1;                      % 1 tune hyper; 0 do not tune
if ~isfield(options,'pseudo')           % number of pseudo inputs
    if options.dataType == 5
        options.pseudo = 49;            % n pseudo inputs in GP
    else
        options.pseudo = 20;            % n pseudo inputs in GP
    end
end
if ~isfield(options,'trainMethod')      % number of training methods for interval censored learning
    options.trainMethod = 1:6;
end
options.evaApprox = 0;          % 1 see how approximation goes
options.wModel = 0;             % 0 GP4C, 1 GP4CW
options.wLow = 1e-6;            % minimum weight in GP4CW
nMethodMax = 6;                 % 1 GP3; 2-4, GP3Int; 5 GP3IntW; 6 LocalEM.
%% initializing the space to store the results
resTest = cell(nMethodMax,options.sampleRepeat);
vBest = inf(nMethodMax,1);      % store the best MISE.
cModel = cell(nMethodMax,1);    % store the best model.

%% repeat the experiment multiple times
for iSampleTimes = 1:options.sampleRepeat 
    options.nSamRepeat = iSampleTimes;
    fprintf('DataType %d\tSampleTimes %d\tPseudo Inputs %d\n',...
        options.dataType,iSampleTimes,options.pseudo);
    % retrieve the data
    data = dataFull{iSampleTimes};
    
    for nMethod = options.trainMethod
        nRun = 1; nCount = 0;    % count the errors
        while nRun == 1          % Safety precautions in case something goes wrong.
            try
                switch nMethod
                    case 1       % complete data with GP
                        options.wModel = 0;
                        [model,nTime] = varTrainGP3Book(data,options);
                    case {2,3,4} % GP4C with different b (2,3,4)=(1,0,0.3)
                        options.wModel = 0;    % GP4C
                        options.inferType = nMethod-1;
                        [model,nTime] = varTrainIntCen(data,options);
                    case 5       % GP4CW
                        options.wModel = 1;    % GP4CW
                        options.inferType = 3;
                        options.nSamRepeat = iSampleTimes;
                        [model,nTime] = varTrainIntCen(data,options);
                    case 6       % LocalEM
                        [model,nTime] = localEM(data,options);
                end
                
                resTmp = testmodel(data,model,options);
                nRun = 0;
            catch ME
                disp(ME.identifier);
                nCount = nCount + 1;
                fprintf('Something weird here');
                if nCount == 3 % wrong three times
                    fprintf('Something weird here');
                    save('dataWrong.mat','data','options');
                    data = generateIntCenBook(options);
                    save('dataReplace.mat','data','options');
                end
            end
       end
        resTmp.nTime = nTime;
        % store the best model
        if resTmp.nMISE < vBest(nMethod)
            vBest(nMethod) = resTmp.nMISE;
            cModel{nMethod} = model;
        end
        displayResult(resTmp);
        resTest{nMethod,iSampleTimes} = resTmp;
    end
    
%     %% censored data with GP
%     for nInfer = arrCan % approximation option
%         nRun = 1;
%         nCount = 0;
%         options.inferType = nInfer;
%         while nRun == 1
%             try
%                 [model,nTime] = varTrainIntCen(data,options);
%                 resTmp = testmodel(data,model,options);
%                 nRun = 0;
%             catch ME
%                 disp(ME.identifier);
%                 nCount = nCount + 1;
%                 fprintf('Something weird here');
%                 if nCount == 3 % wrong three times
%                     fprintf('Something weird here');
%                     save('dataWrong.mat','data','options');
%                     data = generateIntCenBook(options);
%                     save('dataReplace.mat','data','options');
%                 end
%             end
%         end
%         resTmp.nTime = nTime;
%         % store the best model
%         if resTmp.nMISE < vBest(nMethod)
%             vBest(nMethod) = resTmp.nMISE;
%             cModel{nMethod} = model;
%         end
%         displayResult(resTmp);
%         resTest{nMethod,iSampleTimes} = resTmp;
%         nMethod = nMethod+1;
%     end
%     
%     %% censored data with GP and adding additional weight
%     if options.dataWeight == 1
%         options.wModel = 1;
%         options.wLow = 1e-6;
%         options.inferType = 3;
%         options.nSamRepeat = iSampleTimes;
%         [model,nTime] = varTrainIntCen(data,options);
%         resTmp = testmodel(data,model,options);
%         resTmp.nTime = nTime;
%         % store the best model
%         if resTmp.nMISE < vBest(nMethod)
%             vBest(nMethod) = resTmp.nMISE;
%             cModel{nMethod} = model;
%         end
%         displayResult(resTmp);
%         resTest{nMethod,iSampleTimes} = resTmp;
%         nMethod = nMethod+1;
%         options.wModel = 0;
%     end
%     
%     %% Benchmark -- localEM
%     if options.bench == 1
%         % run two sets of accuracy
%         [model,nTime] = localEM(data,options);
%         resTmp = testmodel(data,model,options);
%         resTmp.nTime = nTime;
%         resTmp.kernel = model.bandwidth;
%         % store the best model
%         if resTmp.nMISE < vBest(nMethod)
%             vBest(nMethod) = resTmp.nMISE;
%             cModel{nMethod} = model;
%         end
%         displayResult(resTmp);
%         resTest{nMethod,iSampleTimes} = resTmp;
%     end
%     fprintf('\n');
end
end
