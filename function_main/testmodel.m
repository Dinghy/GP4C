function resTest = testmodel(data,model,options)
%% test the model
% Input :  data      -- struct, required fields
%          model     -- struct, required fields
%               Benchmark: model.pre
%               True     : []
%               GP       : model.var.Sigma,model.var.Mu,model.var.L,model.GP.logtheta
%          options   -- struct, required fields
%               options.testSample
% Output:  resTest   -- struct
%          field:
%               nMISE       -- MISE
%               nSampMISE   -- MISE when using sampling
%               nTestLike   -- test likelihood
%               vTestLike   -- test likelihood for each subject
%               nSampLike   -- test likelihood when using sampling
%               vSampLike   -- test likelihood when using sampling for each subject
%               nTrainLike  -- train likelihood
%               vTrainLike  -- train likelihood for each subject

%% High resolution grid
vT = linspace(data.XBeg,data.XEnd,options.testSample)';
%% Basic computation for variational inference method
% do the basic computation testmodelSave
if isfield(model,'var')
    vLogtheta = model.GP.logtheta;
    vMu = model.var.Mu;
    Kmm = feval(model.cov{:},vLogtheta,model.Xm);
    Knm = feval(model.cov{2}{1},vLogtheta(1:end-1),vT,model.Xm);
    Knn = feval(model.cov{:},vLogtheta,vT,vT);
    KnmInvKmm = Knm/Kmm;
    
    % calculate the integral for each subject in the test data set
    if options.wModel == 1  
        nGamma = exp(2*model.GP.logtheta(length(model.Tend)+1));
        mTmpSigmaMu = Kmm\(model.var.Sigma+model.var.Mu*model.var.Mu')-eye(model.m);
        % on test set
        vfInt = zeros(data.nTest,1);
        for u = 1:data.nTest
            nBeg = min(data.centest{u}.vCInt(:,1));
            nEnd = max(data.centest{u}.vCInt(:,2));
            mPhi = MatrixPhi(model.GP.logtheta,model,nBeg,nEnd);
            invKmmPhi = Kmm\mPhi;
            vfInt(u) = nGamma*prod(nEnd-nBeg)+trace(invKmmPhi*mTmpSigmaMu);
        end
        vEstWTest = cellfun(@(x)sum(x.vCen),data.centest)'./vfInt;
        vEstWTest = max(options.wLow,vEstWTest);
        
        % on train set
        vfInt = zeros(data.nTrain,1);
        for u = 1:data.nTrain
            nBeg = min(data.censor{u}.vCInt(:,1));
            nEnd = max(data.censor{u}.vCInt(:,2));
            mPhi = MatrixPhi(model.GP.logtheta,model,nBeg,nEnd);
            invKmmPhi = Kmm\mPhi;
            vfInt(u) = nGamma*prod(nEnd-nBeg)+trace(invKmmPhi*mTmpSigmaMu);
        end
        vEstWTrain = cellfun(@(x)sum(x.vCen),data.censor)'./vfInt;
        vEstWTrain = max(options.wLow,vEstWTrain);
    end
    
    mSigma = 1e-10*eye(length(vT))+Knn-KnmInvKmm*Knm';
    mL = chol(mSigma);
    mSamplePre = zeros(length(vT),options.testSampTimes);
    for is = 1:options.testSampTimes
        vfM = vMu + model.var.L*randn(size(vMu));
        vfN = KnmInvKmm*vfM+mL*randn(size(vT));
        mSamplePre(:,is) = vfN.^2;
    end
    vPreMean = (KnmInvKmm*vMu).^2+...
        sum(exp(2*vLogtheta(end-1:end)))+...
        sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);
    
    if options.wModel == 1
        f.x = vT;
        mEstWTest = ones(data.nTest,options.testSampTimes);
        mEstWTrain = ones(data.nTrain,options.testSampTimes);
        for is = 1:options.testSampTimes
            f.y = mSamplePre(:,is);
            % on test set
            for u = 1:data.nTest
                nBeg = min(data.centest{u}.vCInt(:,1));
                nEnd = max(data.centest{u}.vCInt(:,2));
                nfInt = SimpsonRule(nBeg,nEnd,[],f);
                mEstWTest(u,is) = max(options.wLow,sum(data.centest{u}.vCen)/nfInt);
            end
            % on train set
            for u = 1:data.nTrain
                nBeg = min(data.censor{u}.vCInt(:,1));
                nEnd = max(data.censor{u}.vCInt(:,2));
                nfInt = SimpsonRule(nBeg,nEnd,[],f);
                mEstWTrain(u,is) = max(options.wLow,sum(data.censor{u}.vCen)/nfInt);
            end
        end
    end
end
%% MISE
nMISE = 0;
if ~isfield(options,'dataWeight')
    options.dataWeight = 0;
end
if isfield(data,'fBasic')
    if isempty(model)                      % true model
        nMISE = 0;
    elseif isfield(model,'pre')            % Benchmark
        if ~iscolumn(model.pre)
            model.pre = model.pre';
        end
        f.x = vT;
        nMISE = 0;
        for u=1:data.nTest
            f.y = (data.centest{u}.fCurrent(vT)-model.pre).^2;
            nMISE = nMISE + SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
        end
        nMISE = nMISE/data.nTest;
    elseif isfield(model,'var')        % Variational inference
        % result by the mean function
        f.x = vT;
        nMISE = 0;
        for u = 1:data.nTest
            % complex model
            if options.wModel == 1
                f.y = (data.centest{u}.fCurrent(vT)-vEstWTest(u)*vPreMean).^2;
            else
                f.y = (data.centest{u}.fCurrent(vT)-vPreMean).^2;
            end
            nMISE = nMISE + SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
        end
        nMISE = nMISE/data.nTest;
        % result by sampling
        f.x = vT;
        nSampMISE = 0;
        for u = 1:data.nTest
            for is = 1:options.testSampTimes
                if options.wModel == 1
                    f.y = (data.centest{u}.fCurrent(vT)-mEstWTest(u,is)*mSamplePre(:,is)).^2;
                else
                    f.y = (data.centest{u}.fCurrent(vT)-mSamplePre(:,is)).^2;
                end
                
                nSampMISE = nSampMISE + SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
            end
        end
        nSampMISE = nSampMISE/options.testSampTimes/data.nTest;
        resTest.nSampMISE = nSampMISE;
    else
        warning('Illegal model in testmodel');
    end
end
resTest.nMISE = nMISE;
%% Test Likelihood
nTestLike = 0;
if isempty(model)                      % true test likelihood
    [nTestLike,vTestLike,nTrainLike,vTrainLike] = calTestLikelihood(data,[],[],[]);
elseif isfield(model,'pre')            % Benchmark
    f.x = vT;
    f.y = model.pre;
    [nTestLike,vTestLike,nTrainLike,vTrainLike] = calTestLikelihood(data,f,[],[]);
elseif isfield(model,'var')            % Variational inference
    % by the result for the mean function
    f.x = vT;
    f.y = vPreMean;
    if options.wModel == 1
        [nTestLike,vTestLike,nTrainLike,vTrainLike] = calTestLikelihood(data,f,vEstWTest,vEstWTrain);
    else
        [nTestLike,vTestLike,nTrainLike,vTrainLike] = calTestLikelihood(data,f,[],[]);
    end
    % by sampling
    f.x = vT;
    mTmp = zeros(length(vTestLike),options.testSampTimes);
    mTmp2 = zeros(length(vTrainLike),options.testSampTimes);
    for is = 1:options.testSampTimes
        f.y = mSamplePre(:,is);
        if options.wModel == 1   % GP4CW
            [~,vTmp,~,vTmp2] = calTestLikelihood(data,f,mEstWTest(:,is),mEstWTrain(:,is)); 
        else                     % GP4C
            [~,vTmp,~,vTmp2] = calTestLikelihood(data,f,[],[]); 
        end
        mTmp(:,is) = vTmp;
        mTmp2(:,is) = vTmp2;
    end
    vMaxTmp = max(mTmp,[],2);
    vSampLike = log(mean(exp(mTmp-repmat(vMaxTmp,1,options.testSampTimes)),2))+vMaxTmp;
    % problematic one here % nSampLike = sum(vSampLike);
    vTmp = sum(mTmp,1);
    nMax = max(vTmp);
    nSampLike = log(mean(exp(vTmp-nMax)))+nMax;
    
    resTest.nSampLike = nSampLike;
    resTest.vSampLike = vSampLike;
    
    vMaxTmp = max(mTmp2,[],2);
    vSampLike = log(mean(exp(mTmp2-repmat(vMaxTmp,1,options.testSampTimes)),2))+vMaxTmp;
    resTest.vSampTrainLike = vSampLike;
else
    warning('Illegal model in testmodel');
end
resTest.nTestLike = nTestLike;
resTest.vTestLike = vTestLike;
resTest.nTrainLike = nTrainLike;
resTest.vTrainLike = vTrainLike;
end

%     if options.dataWeight == 0                 % normal data with single basis intensity
%         if isempty(model)                      % true model
%             nMISE = 0;
%         elseif isfield(model,'pre')            % Benchmark
%             if ~iscolumn(model.pre)
%                 model.pre = model.pre';
%             end
%             f.x = vT;
%             f.y = (data.fBasic(vT)-model.pre).^2;
%             nMISE = SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
%         elseif isfield(model,'var')        % Variational inference
%             % result by the mean function
%             f.x = vT;
%             f.y = (data.fBasic(vT)-vPreMean).^2;
%             nMISE = SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
%             % result by sampling
%             f.x = vT;
%             nSampMISE = 0;
%             for is = 1:options.testSampTimes
%                 f.y = (data.fBasic(vT)-mSamplePre(:,is)).^2;
%                 nSampMISE = nSampMISE +SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
%             end
%             nSampMISE = nSampMISE/options.testSampTimes;
%             resTest.nSampMISE = nSampMISE;
%         else
%             warning('Illegal model in testmodel');
%         end
%     else               % multiple basis intensity function