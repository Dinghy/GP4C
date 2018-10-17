function [nMISE,nTestLike] = testMSE(data,model,options)
% Required Field:
%   -- options.testSample
%   -- options.testSampTimes

% rewrite the function for evaluating the intensity function
nMISE = 0;      % MISE
nTestLike = 0;  % Test Likelihood

vT = linspace(data.XBeg,data.XEnd,options.testSample)';
if isfield(data,'fBasic')
    if isempty(model)                      % true model
        for u = 1:data.nTest
            nTestLike = nTestLike + sum(log(data.fBasic(data.Test{u}.X)))-...
                SimpsonRule(data.XBeg,data.XEnd,options.testSample,data.fBasic);
        end
    elseif isfield(model,'flam')            % Benchmark
        f.x = vT;
        for u=1:data.nTest
            f.y = (data.fBasic(vT)-model.flam(vT)).^2;
            nMISE = nMISE + SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
            f.y = model.flam(vT);
            nTestLike = nTestLike + sum(log(model.flam(data.test{u}.X)))-...
                SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
        end
        nMISE = nMISE/data.nTest;
    elseif isfield(model,'var')        % Variational inference
        % result by the mean function
        nMISE = 0;
        for u = 1:data.nTest
            f.x = vT;
            [vMeanPre,mSamplePre] = gpPredict(model,[vT;data.test{u}.X],options);
            % complex model
            f.y = (data.fBasic(vT)-vMeanPre(1:length(vT))).^2;
            nMISE = nMISE + SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
            vtest = zeros(1,options.testSampTimes);
            for is = 1:options.testSampTimes
                f.y = mSamplePre(1:length(vT),is);
                nInt = SimpsonRule(data.XBeg,data.XEnd,options.testSample,f);
                vtest(is) = sum(log(mSamplePre(length(vT)+1:end,is)))-nInt;
            end
            vtest_exp = exp(vtest-max(vtest));
            nTestLike = log(mean(vtest_exp))+max(vtest);
        end
        nMISE = nMISE/data.nTest;
    else
        warning('Illegal model in testmodel');
    end
end
end

function [vMeanPre,mSamplePre] = gpPredict(model,vT,options)
% Perform the prediction by the GP method
% Input:
%   model: model that contains the parameters
%   vN   : prediction points
% Output:
%   vMeanPre    : Prediction of the mean intensity function
%   mSamplePre  : Prediction of the sampled intensity functions
vLogtheta = model.GP.logtheta;

vMu = model.var.Mu;
Kmm = feval(model.cov{:},vLogtheta,model.Xm);
Knm = feval(model.cov{2}{1},vLogtheta(1:end-1),vT,model.Xm);
Knn = feval(model.cov{:},vLogtheta,vT);
KnmInvKmm = Knm/Kmm;

vMeanPre = (KnmInvKmm*vMu).^2+...
    sum(exp(2*vLogtheta(end-1:end)))+...
    sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);

mSigma = 1e-10*eye(length(vT))+Knn-KnmInvKmm*Knm';
mL = chol(mSigma);
mSamplePre = zeros(length(vT),options.testSampTimes);
for is = 1:options.testSampTimes
    vfM = vMu + model.var.L*randn(size(vMu));
    vfN = KnmInvKmm*vfM+mL*randn(size(vT));
    mSamplePre(:,is) = vfN.^2;
end
end