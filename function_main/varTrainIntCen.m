function [model,nTime] = varTrainIntCen(data,options)
% Difference from varTrainIntCen
%       1. Assume all experiments come from the same intensity function
%       2. Assume all experiments from the same module but with different
%          amplitude

%% Check the field of the training
if ~isfield(options,'optOffset')  % whether to optimize the offset
    options.optOffset = 1;        % default
end
if ~isfield(options,'detailPlot') % whether to do the plots
    options.detailPlot = 0;
end

if ~isfield(options,'offset')     % offset in approximation
    options.offset = 0.30;
    switch options.inferType 
        case 1
            options.offset = 1;
        case 2
            options.offset = 0;
        case 3
            options.offset = 0.3;
        case 4
            % options.inferType = 3;
            options.offset = 0.7;
    end
end
if options.figure == 1
    fig = figure;
    set(fig,'Position',[0,0,700,500]);
end
if options.detailPlot == 1
    fig2 = figure;
    para.Mu = zeros(1,options.maxIter);
    para.L = zeros(1,options.maxIter);
    para.theta = zeros(length(data.XBeg)+2,options.maxIter);
end

%% initialization
model = varInit(data,options);
load('gzdgz.mat');
model.gz=gz;
model.dgz=dgz;
%% main iteration
nTime = 0;
model.cenEva = [];
[LB,UB] = optConstraint(model);
tStart = tic;
nELBOpre = -Inf;
iIter = 1;
% only save a local copy of current u
modelD = preKernel(data,model);
while iIter <= options.maxIter
    if options.detailPlot == 1
        para.Mu(iIter) = sum(abs(model.var.Mu));
        para.L(iIter) = sum(sum(abs(model.var.L)));
        para.theta(:,iIter) = model.GP.logtheta;
    end
    modelD = preVar(data,model,modelD);
    if options.evaApprox == 1
        model = evaluateApprox(data,model,modelD);
    end
    if options.debug == 1
        V = extractMuSigma(model);
        diff = checkGrad('calDVar',V,data,model,modelD);
        disp(diff);
        V = extractTheta(model,options);
        diff = checkGrad('calDKer',V,data,model,modelD,options);
        disp(diff);
        error('Debug Finished');
    end
    %% optimize mu,sigma
    vVar=extractMuSigma(model);
    myObj = @(VV)calDVar(VV,data,model,modelD);
    [V,nELBO] = minConf_TMP(myObj,vVar,LB,UB,struct('verbose',0,'maxIter',10));
    % return the parameter
    model = returnMuSigma(model,V);
    modelD = preVar(data,model,modelD); % update
    %% update the weight parameter if needed
    if options.wModel == 1
        model.var.W = cellfun(@(x)sum(x),model.st.count)./modelD.vfInt;
        model.var.W = max(model.wLow,model.var.W);  % threshold at the lowest value
    end
    %% optimize hyper
    if options.hyper == 1
        vTheta = extractTheta(model,options);
        myObj = @(VV)calDKer(VV,data,model,modelD,options);
        [V,~] = minConf_TMP(myObj,vTheta,...
            model.GP.cons(1,1:length(vTheta))',model.GP.cons(2,1:length(vTheta))',struct('verbose',0,'maxIter',10));
        model = returnTheta(model,V,options);
        modelD = preKernel(data,model);
        modelD = preVar(data,model,modelD); % update
    end
    %% Plotting part
    if options.figure == 1
        if options.evaApprox == 1
            fprintf('Correlation 1:%.4f\t 2:%.4f\t 3:%.4f\n',cov(model.cenEva(2,:)-model.cenEva(1,:)),...
                cov(model.cenEva(3,:)-model.cenEva(1,:)),cov(model.cenEva(4,:)-model.cenEva(1,:)));
            plotResult(model,modelD,options,data,fig,1);
        else
            plotResult(model,modelD,options,data,fig,2);
        end
    end
    %% break if finished
    if abs(nELBO-nELBOpre)<=1e-4*abs(nELBO)
        nTime = toc(tStart);
        break
    else
        % fprintf('Iteration %d\t:%.4f\tTrain:%.4f\tTest:%.4f\n',iIter,nELBO,nTrain,nTest);
        nELBOpre = nELBO;
        iIter = iIter + 1;
    end
end
if options.detailPlot == 1
    plotDetail(fig2,para,iIter);
end
end

function model = varInit(data,options)
%% Initialization part for the model
model.const.C = 0.5772156649;          % Euler's constant
model.m = options.pseudo;                          % number of pseudo inputs
model.Tstart = data.XBeg;              % time start
model.Tend = data.XEnd;                % time end
model.cov = {'covSum',{'covSEard','covNoise'}};


model.Xm = linspace(model.Tstart,model.Tend,model.m)';
model.plotXm = linspace(model.Tstart,model.Tend,100)';

model.GP.cons = zeros(1,2);
% model.inferType = options.inferType;
model.offset = options.offset;
model.wModel = options.wModel;

% initialization and counting
% space allocation
model.NumData = zeros(data.nU,1);
model.st.period = cell(data.nU,1);
model.st.range = zeros(data.nU,2);
model.st.count = cell(data.nU,1);
for u=1:data.nU
    model.NumData(u) = length(data.censor{u}.X)+sum(data.censor{u}.vCen);  % number of points in each dataset
    model.st.range(u,:) = [data.censor{u}.vCInt(1,1),data.censor{u}.vCInt(end,2)];
    btag = data.censor{u}.vCen>0;
    model.st.period{u} = data.censor{u}.vCInt(btag,:);
    model.st.count{u} = data.censor{u}.vCen(btag);
end
% model.prior.gbase = sqrt(mean(model.NumData)/(model.Tend-model.Tstart));
model.prior.gbase = sqrt(mean(model.NumData./(model.st.range(:,2)-model.st.range(:,1))));
% hyper-parameter for GP
D = length(data.XBeg);
if options.hyper == 0
    model.GP.logtheta = [1.5136;-0.4551;-6.9078];%[log(2);0.5*log(0.1);log(1e-3)];
else
    dd = log((model.Tend-model.Tstart)')/2;
    model.GP.logtheta = zeros(D+2,1);
    model.GP.logtheta(1:D) = dd;
    model.GP.logtheta(D+1) = 0.5*log(0.5*model.prior.gbase^2);
    model.GP.logtheta(D+2) = log(1e-3);
end
model.GP.cons(1,1:D) = -inf;
model.GP.cons(2,1:D) = inf;
model.GP.cons(1,D+1) = -inf;
model.GP.cons(2,D+1) = inf;
model.GP.cons(1,D+2) = log(1e-2);
model.GP.cons(2,D+2) = log(1e-5);

% Variational distribution
if options.wModel == 1
    model.wLow = options.wLow;
    vCount = cellfun(@(x)sum(x),model.st.count);
    vTime = model.st.range(:,2)-model.st.range(:,1);
    model.var.W = vCount./vTime;
    model.var.W = max(model.wLow,model.var.W);
end

model.var.Mu = sqrt(0.5*model.prior.gbase^2)*(0.1+rand(model.m,1)); 


model.var.Sigma = feval(model.cov{:},model.GP.logtheta,model.Xm);
model.var.L = chol(model.var.Sigma,'lower');
end

function vV = extractTheta(model,options)
%% extract theta
switch options.hyperOpt
    case 0 % do not optimize the jitter part
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta(1:end-1);
        else
            vV = model.GP.logtheta(1:end-1)';
        end
    case 1 % do not optimize the length term in a GP
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta(1:end-2);
        else
            vV = model.GP.logtheta(1:end-2)';
        end
    case 2  % optimize all hyper-parameters in a GP
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta;
        else
            vV = model.GP.logtheta';
        end
end
end

function model = returnTheta(model,V,options)
% return theta
switch options.hyperOpt
    case 0
        N = length(model.GP.logtheta)-1;
        model.GP.logtheta(1:end-1) = V(1:N);
    case 1
        N = length(model.GP.logtheta)-2;
        model.GP.logtheta(1:end-2) = V(1:N);
    case 2
        model.GP.logtheta= V;
end
end

function vVar = extractMuSigma(model)
%% extract Mu, Sigma
vVar=[model.var.Mu;model.var.L(tril(true(size(model.var.L))))];
end

function model = returnMuSigma(model,V)
%% return Mu, Sigma
M = model.m;
M2= M*(M+1)/2;
nStart = 1;
model.var.Mu = V(nStart:nStart+M-1);
mtmp = tril(ones(M));
nStart = nStart+M;
mtmp(mtmp==1) = V(nStart:nStart+M2-1);
model.var.L = mtmp;
model.var.Sigma = mtmp*mtmp';
end

function modelD = preKernel(data,model,modelD)
%  PRE-CALCULATION FOR Parameters that are only related to kernel
%  hyper-parameters, Kmm, mPhi, invKmmPhi, diagKnn
if nargin == 2 % initialization
    modelD.numPar = length(model.GP.logtheta);
end
modelD.diagKnn = sum(exp(2*model.GP.logtheta(end-1:end)));
modelD.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xm);
modelD.KmmDx = cell(1,modelD.numPar);
for k=1:modelD.numPar
    modelD.KmmDx{k} = feval(model.cov{:}, model.GP.logtheta, model.Xm, [], k);
end

% calculate phi for each censored interval
if ~isfield(modelD,'DmCenPhi') % initialization
    % entire integral
    modelD.invKmmPhi = cell(data.nU,1);% modelD.Kmm\mPhi;
    modelD.DmPhidx = cell(data.nU,1); % DmPhidx;
    % censored integral
    modelD.DmCenPhi = cell(data.nU,1);
    modelD.mCenInvKmmPhi = cell(data.nU,1);
    for u=1:data.nU
        nCen = length(model.st.count{u});
        modelD.DmCenPhi{u} = cell(nCen,1);
        modelD.mCenInvKmmPhi{u} = cell(nCen,1);
    end
end

mapObj = containers.Map('KeyType','double','ValueType','any');
for u = 1:data.nU
    % entire integral
    nBeg = model.st.range(u,1);
    nEnd = model.st.range(u,2);
    % look up in the map
    ntmp = nBeg*(data.XEnd+1)+nEnd;
    if isKey(mapObj,ntmp)
        c = mapObj(ntmp);
        invKmmPhi = c.invKmmPhi;
        DmPhidx = c.DmPhidx;
    else
        [mPhi,DmPhidx] = MatrixPhi(model.GP.logtheta,model,nBeg,nEnd);
        invKmmPhi = modelD.Kmm\mPhi;
        c.invKmmPhi = invKmmPhi; 
        c.DmPhidx = DmPhidx;
        mapObj(ntmp) = c;
    end
    modelD.DmPhidx{u} = DmPhidx;
    modelD.invKmmPhi{u} = invKmmPhi;
    % censored integral
    for i=1:length(model.st.count{u})
        nBeg = model.st.period{u}(i,1);
        nEnd = model.st.period{u}(i,2);
        % look up in the map
        ntmp = nBeg*(data.XEnd+1)+nEnd;
        if isKey(mapObj,ntmp)
            c = mapObj(ntmp);
            invKmmPhi = c.invKmmPhi;
            DmPhidx = c.DmPhidx;
        else
            [mPhi,DmPhidx] = MatrixPhi(model.GP.logtheta,model,nBeg,nEnd);
            invKmmPhi = modelD.Kmm\mPhi;
            c.invKmmPhi = invKmmPhi;
            c.DmPhidx = DmPhidx;
            mapObj(ntmp) = c;
        end
        modelD.DmCenPhi{u}{i} = DmPhidx;
        modelD.mCenInvKmmPhi{u}{i} = invKmmPhi;
    end
end
end

function modelD = preVar(data,model,modelD)
%  PRE-CALCULATION FOR Parameters used in update for m,S
D=length(model.Tend);
% common parts
modelD.InvKmmMu = modelD.Kmm\model.var.Mu;
modelD.InvKmmSigma = modelD.Kmm\model.var.Sigma;
modelD.InvKmmL = modelD.Kmm\model.var.L;
nGamma = exp(2*model.GP.logtheta(D+1));
InvKmmSigmaMu = modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu'-eye(model.m);
% integral for the whole space
modelD.vfInt = zeros(data.nU,1);

for u = 1:data.nU
    modelD.vfInt(u) = nGamma*prod(model.st.range(u,2)-model.st.range(u,1))+trace(modelD.invKmmPhi{u}*InvKmmSigmaMu);
end
% individual user
if ~isfield(modelD,'vfIntCen') % initialization
    modelD.vfIntCen = cell(data.nU,1);
    for u=1:data.nU
        modelD.vfIntCen{u} = ones(length(model.st.count{u}),1);
    end
end

% precomputation of the approximation of integral
InvKmmMu = modelD.InvKmmMu*model.var.Mu';
InvKmmMuSig = InvKmmMu+model.offset*modelD.InvKmmSigma-...
    model.offset*eye(model.m);

for u=1:data.nU
    nCen = length(model.st.count{u});
    %     switch model.inferType
    %         case 1
    %             for i=1:nCen
    %                 nBeg = model.st.period{u}(i,1);
    %                 nEnd = model.st.period{u}(i,2);
    %                 modelD.vfIntCen{u}(i) = nGamma*prod(nEnd-nBeg)+trace(modelD.mCenInvKmmPhi{u}{i}*InvKmmSigmaMu);
    %             end
    %         case 2
    %             for i=1:nCen
    %                 modelD.vfIntCen{u}(i) = trace(modelD.mCenInvKmmPhi{u}{i}*InvKmmMu);
    %             end
    %         case 3
    for i=1:nCen
        nBeg = model.st.period{u}(i,1);
        nEnd = model.st.period{u}(i,2);
        modelD.vfIntCen{u}(i) = trace(modelD.mCenInvKmmPhi{u}{i}*InvKmmMuSig)+...
            model.offset*nGamma*prod(nEnd-nBeg);
    end
    %     end
end
end

function [f,df] = calDVar(V,data,model,modelD)
% calculate the derivative and the current likelihood
model = returnMuSigma(model,V);
modelD = preVar(data,model,modelD);
modelD.f = 0;
modelD.var.Mu = zeros(size(model.var.Mu));
modelD.var.L = zeros(size(model.var.L));
%% add irrelevant part (weight part if wModel = 1, integral part)
if model.wModel == 1 % GP4CW 
    for u = 1:data.nU
        modelD.f = modelD.f + sum(model.st.count{u})*log(model.var.W(u));
        modelD.f = modelD.f-model.var.W(u)*modelD.vfInt(u);
        modelD.var.Mu = modelD.var.Mu-2*model.var.W(u)*modelD.invKmmPhi{u}*modelD.InvKmmMu;
        modelD.var.L = modelD.var.L-2*model.var.W(u)*modelD.invKmmPhi{u}*modelD.InvKmmL;
    end
else                 % GP4C
    for u = 1:data.nU
        modelD.f = modelD.f-modelD.vfInt(u);
        modelD.var.Mu = modelD.var.Mu-2*modelD.invKmmPhi{u}*modelD.InvKmmMu;
        modelD.var.L = modelD.var.L-2*modelD.invKmmPhi{u}*modelD.InvKmmL;
    end
end
%% add the interval-censored part
mInvKmmPhi = zeros(model.m);
for u = 1:data.nU
    if ~isempty(model.st.count{u})
        modelD.f = modelD.f + sum(model.st.count{u}.*log(modelD.vfIntCen{u}));
        for ic = 1:length(model.st.count{u})
            mInvKmmPhi = mInvKmmPhi + model.st.count{u}(ic)/modelD.vfIntCen{u}(ic)*modelD.mCenInvKmmPhi{u}{ic};
        end
    end
end
% switch model.inferType
%     case 1 % upper bound
%         modelD.var.Mu = modelD.var.Mu + 2*mInvKmmPhi*modelD.InvKmmMu;
%         modelD.var.L = modelD.var.L + 2*mInvKmmPhi*modelD.InvKmmL;
%     case 2 % lower bound
%         modelD.var.Mu = modelD.var.Mu + 2*mInvKmmPhi*modelD.InvKmmMu;
%     case 3 % combination
modelD.var.Mu = modelD.var.Mu + 2*mInvKmmPhi*modelD.InvKmmMu;
modelD.var.L = modelD.var.L + 2*model.offset*mInvKmmPhi*modelD.InvKmmL;
% end
%% KL divergence and gradient
InvKmmMuG = modelD.Kmm\(model.var.Mu-model.prior.gbase);
InvKmmSigmaMuG=modelD.InvKmmSigma+InvKmmMuG*(model.var.Mu-model.prior.gbase)';
modelD.f = modelD.f+0.5*(model.m-trace(InvKmmSigmaMuG))...
    +sum(log(diag(model.var.L)))-0.5*logdet(modelD.Kmm);
modelD.var.Mu = modelD.var.Mu-InvKmmMuG;
modelD.var.L = modelD.var.L+(diag(1./diag(model.var.L))-modelD.InvKmmL);
%% extract the gradient
f = - modelD.f;
df = - extractMuSigma(modelD);
end

function [f,df] = calDKer(W,data,model,modelD,options)
% calculate the function value and gradient
% return the parameters in W
model = returnTheta(model,W,options);

modelD = preKernel(data,model,modelD);
modelD = preVar(data,model,modelD);
%% initializing and calculation
modelD.f=0;
modelD.GP.logtheta=zeros(size(model.GP.logtheta));
modelD.InvKmmSigmaInvKmm = modelD.InvKmmSigma/modelD.Kmm;
%% Integration over Tlim
nGamma = exp(2*model.GP.logtheta(length(model.Tend)+1));
InvKmmSigmaMu=modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu';
mtmp_a = (eye(model.m)-InvKmmSigmaMu)/modelD.Kmm;
mtmp_b = zeros(model.m);
if options.wModel == 1 % GP4CW
    for u = 1:data.nU
        modelD.f = modelD.f-model.var.W(u)*modelD.vfInt(u);
        mtmp_b = mtmp_b + model.var.W(u)*(2*modelD.invKmmPhi{u}*InvKmmSigmaMu-modelD.invKmmPhi{u});
        modelD.GP.logtheta = modelD.GP.logtheta+...
            model.var.W(u)*cellfun(@(x)trace(mtmp_a*x),modelD.DmPhidx{u})';
    end
    if length(W) > length(model.Tend)
        modelD.GP.logtheta(end-1) = modelD.GP.logtheta(end-1)-...
            nGamma*2*sum(model.var.W.*(model.st.range(:,2)-model.st.range(:,1)));
    end
else                   % GP4C
    for u = 1:data.nU
        modelD.f = modelD.f-modelD.vfInt(u);
        mtmp_b = mtmp_b + (2*modelD.invKmmPhi{u}*InvKmmSigmaMu-modelD.invKmmPhi{u});
        modelD.GP.logtheta = modelD.GP.logtheta+...
            cellfun(@(x)trace(mtmp_a*x),modelD.DmPhidx{u})';
    end
    if length(W) > length(model.Tend)
        modelD.GP.logtheta(end-1) = modelD.GP.logtheta(end-1)-nGamma*2*sum(model.st.range(:,2)-model.st.range(:,1));
    end
end

%% add the censored part
% switch model.inferType
%     case 1
%         mtmp_a = (eye(model.m)-InvKmmSigmaMu)/modelD.Kmm;
%     case 2
%         InvKmmMu=modelD.InvKmmMu*model.var.Mu';
%         mtmp_a = -InvKmmMu/modelD.Kmm;
%     case 3
InvKmmbSigmaMu = model.offset*modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu';
mtmp_a = (model.offset*eye(model.m)-InvKmmbSigmaMu)/modelD.Kmm;
% end
for u = 1:data.nU
    if ~isempty(model.st.count{u})
        modelD.f = modelD.f + sum(model.st.count{u}.*log(modelD.vfIntCen{u}));
        nCen = length(model.st.count{u});
        %         switch model.inferType
        %             case 1 % upper bound
        %                 for i = 1:nCen
        %                     mtmp_b = mtmp_b - ...
        %                         model.st.count{u}(i)/modelD.vfIntCen{u}(i)*(2*modelD.mCenInvKmmPhi{u}{i}*InvKmmSigmaMu-modelD.mCenInvKmmPhi{u}{i});
        %                     modelD.GP.logtheta = modelD.GP.logtheta-model.st.count{u}(i)/modelD.vfIntCen{u}(i)*...
        %                         cellfun(@(x)trace(mtmp_a*x),modelD.DmCenPhi{u}{i})';
        %                     if length(W) > length(model.Tend)
        %                         nT = model.st.period{u}(i,2)-model.st.period{u}(i,1);
        %                         modelD.GP.logtheta(end-1) = modelD.GP.logtheta(end-1)+model.st.count{u}(i)/modelD.vfIntCen{u}(i)*nGamma*2*nT;
        %                     end
        %                 end
        %             case 2 % lower bound
        %                 for i = 1:nCen
        %                     mtmp_b = mtmp_b - ...
        %                         model.st.count{u}(i)/modelD.vfIntCen{u}(i)*2*modelD.mCenInvKmmPhi{u}{i}*InvKmmMu;
        %                     modelD.GP.logtheta = modelD.GP.logtheta-model.st.count{u}(i)/modelD.vfIntCen{u}(i)*...
        %                         cellfun(@(x)trace(mtmp_a*x),modelD.DmCenPhi{u}{i})';
        %                 end
        %             case 3 % combination
        for i = 1:nCen
            mtmp_b = mtmp_b - ...
                model.st.count{u}(i)/modelD.vfIntCen{u}(i)*(2*modelD.mCenInvKmmPhi{u}{i}*InvKmmbSigmaMu-model.offset*modelD.mCenInvKmmPhi{u}{i});
            modelD.GP.logtheta = modelD.GP.logtheta-model.st.count{u}(i)/modelD.vfIntCen{u}(i)*...
                cellfun(@(x)trace(mtmp_a*x),modelD.DmCenPhi{u}{i})';
            if length(W) > length(model.Tend)
                nT = model.st.period{u}(i,2)-model.st.period{u}(i,1);
                modelD.GP.logtheta(end-1) = modelD.GP.logtheta(end-1)+model.offset*model.st.count{u}(i)/modelD.vfIntCen{u}(i)*nGamma*2*nT;
            end
        end
        %         end
    end
end
mtmp_b = mtmp_b/modelD.Kmm;
modelD.GP.logtheta = modelD.GP.logtheta + cellfun(@(x)trace(mtmp_b*x),modelD.KmmDx)';
%% KL divergence and gradient
InvKmmMuG = modelD.Kmm\(model.var.Mu-model.prior.gbase);
InvKmmSigmaMuG=modelD.InvKmmSigma+InvKmmMuG*(model.var.Mu-model.prior.gbase)';

modelD.f = modelD.f+0.5*(model.m-trace(InvKmmSigmaMuG))...
    +sum(log(diag(model.var.L)))-0.5*logdet(modelD.Kmm);
mtmp = InvKmmSigmaMuG/modelD.Kmm;
modelD.GP.logtheta = modelD.GP.logtheta-...
    0.5*cellfun(@(x)trace(modelD.Kmm\x-mtmp*x),modelD.KmmDx)';
%% extract the variable
f = - modelD.f;
df = - extractTheta(modelD,options);
end

function plotResult(model,modelD,options,data,fig,bchoice)
% plot latent function and confidence interval
ilinewidth = 1.4;
figure(fig);set(gca,'FontSize',13);

switch bchoice
    case 1
        nxmin = min(model.cenEva(:,1));
        nxmax = max(model.cenEva(:,1));
        vx = linspace(nxmin,nxmax,100);
        loglog(vx,vx,'g','LineWidth',ilinewidth);hold on;
        loglog(model.cenEva(:,1),model.cenEva(:,2),'b*','LineWidth',ilinewidth);
        loglog(model.cenEva(:,1),model.cenEva(:,3),'r*','LineWidth',ilinewidth);
        loglog(model.cenEva(:,1),model.cenEva(:,4),'m*','LineWidth',ilinewidth);hold off;
        xlabel('True Integral value by Simpson rule');
        ylabel('Approximation value');
        legend('Ideal','ln(\mu^2+\sigma^2)','ln(\mu^2)','ln(\mu^2+b\sigma^2)','Location','best');
    case 2
        vT = model.plotXm;
        Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
        
        KnmInvKmm = Knm/modelD.Kmm;
        vMean = KnmInvKmm*model.var.Mu;
        vVar = modelD.diagKnn+sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma-eye(model.m))),2);

        vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
        vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
        vMean = vMean.^2+vVar;
        figure(fig);
        subplot('Position',[0.1,0.12,0.8,0.8]);
        if isfield(data,'fBasic') && options.dataWeight == 0
            plot(vT,data.fBasic(vT),'k','LineWidth',ilinewidth);hold on;
            plot(vT,vMean,'r','LineWidth',ilinewidth);
        else
            plot(vT,vMean,'r','LineWidth',ilinewidth);hold on;
        end
        plot(vT,vHigh,'r--','LineWidth',ilinewidth);
        plot(vT,vLow,'r--','LineWidth',ilinewidth);hold off;
        xlim([data.XBeg,data.XEnd]);
        set(gca,'XTick',[]);
        set(gca,'FontSize',12);
end
drawnow;
end

function plotDetail(fig,para,nIter)
% Plot the details of the change of parameters
figure(fig);
subplot(2,3,1);plot(1:nIter,para.Mu(1:nIter));title('Mu');
subplot(2,3,2);plot(1:nIter,para.L(1:nIter));title('L');
subplot(2,3,3);plot(1:nIter,para.theta(1,1:nIter));title('Theta 1');
subplot(2,3,4);plot(1:nIter,para.theta(2,1:nIter));title('Theta 2');
subplot(2,3,5);plot(1:nIter,para.theta(3,1:nIter));title('Theta 3');
end

function norm = checkGrad(objFunc,vInput,varargin)
% DEBUG: function to check whether gradient caculation is correct
%
nEps=1e-3;
norm = 0;
[f,df]=feval(objFunc,vInput,varargin{:});

for i=1:length(vInput)
    vInput(i) = vInput(i) - 0.5*nEps;
    E0 = feval(objFunc,vInput,varargin{:});
    vInput(i) = vInput(i) + nEps;
    E1 =feval(objFunc,vInput,varargin{:});
    vInput(i) = vInput(i) - 0.5*nEps;
    diff_grad = (E1 - E0)/nEps;
    norm = norm + abs(diff_grad - df(i));
    fprintf('%d: %f vs %f: Diff %f\n',i, diff_grad, df(i),abs(diff_grad - df(i)));
end
fprintf('Likelihood: %f\n',f);
end

function model = evaluateApprox(data,model,modelD)
%% for each censored period, calculate the expectation of the integral and the approximated value
np = 30;
ns = 300;
D=length(model.Tend);
nGamma = exp(2*model.GP.logtheta(D+1));
InvKmmMu=modelD.InvKmmMu*model.var.Mu';
InvKmmSigmaMu=modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu';

for u = 1:data.nU
    for ic = 1:data.censor{u}.nCInt
        nBeg = data.censor{u}.vCInt(ic,1);
        nEnd = data.censor{u}.vCInt(ic,2);
        % approximate by simpson's rule
        vx = linspace(nBeg,nEnd,np)';
        nh = vx(2)-vx(1);
        nInt = 0;
        nApp1 = nGamma*prod(nEnd-nBeg)-trace(modelD.mCenInvKmmPhi{u}{ic})+...
            trace(modelD.mCenInvKmmPhi{u}{ic}*InvKmmSigmaMu); % mu^2+sig^2
        nApp2 = trace(modelD.mCenInvKmmPhi{u}{ic}*InvKmmMu);  % mu^2
        nApp3 = trace(modelD.mCenInvKmmPhi{u}{ic}*InvKmmMu)+...
            model.offset*(nGamma*prod(nEnd-nBeg)-trace(modelD.mCenInvKmmPhi{u}{ic})+...
            trace(modelD.mCenInvKmmPhi{u}{ic}*modelD.InvKmmSigma)); % mu^2+b*sig^2
        Knm = feval(model.cov{:},model.GP.logtheta,vx,model.Xm);
        KnmInvKmm = Knm/modelD.Kmm;
        vSigma = sum(exp(2*model.GP.logtheta(end-1:end)))-sum(KnmInvKmm.*Knm,2);
        vSigma = sqrt(vSigma);
        for is = 1:ns
            vfM = model.var.Mu;% + mL*randn(size(vMu));
            vf = KnmInvKmm*vfM+vSigma.*randn(size(vx));
            nInt = nInt + nh/3*(vf(1)^2+vf(np)^2+2*sum(vf(2:2:np-2).^2)...
                +4*sum(vf(1:2:np-1).^2));
        end
        nInt = nInt/ns;
        model.cenEva = [model.cenEva;abs(nInt),abs(nApp1),abs(nApp2),abs(nApp3)];
    end
end
end

