function [model,nTime] = localEM(data,options)
% run the localEM algorithm based on R code
% Input   data    -- basic data structure
%         options -- 
% Output  model   -- return the model, contain model.pre stores the
%                    prediction
%         nTime   -- training time

%% Basic settings and checkings
if ~isfield(options,'nCross')
    options.nCross = 5;
end
if ~isfield(options,'benchPos')
    options.benchPos = 15;
end
vBwSeq = linspace(0.02,sqrt(data.XEnd-data.XBeg),options.benchPos);
%% Set a timer for the training process
tStart = tic;
% split the training data into five parts
nPart = floor(data.nU/options.nCross);
cPart = cell(1,options.nCross);
for i = 1:options.nCross
    cPart{i} = (i-1)*nPart+1:i*nPart;
    if i == options.nCross
        cPart{i} = (i-1)*nPart+1:data.nU;
    end
end
vLikeBw = zeros(size(vBwSeq));
vSaveTest = zeros(data.nU,1);
fprintf('Benchmark progress: ');
for i = 1:length(vBwSeq)
    fprintf('%d\t',i);
    modelD.bandwidth = vBwSeq(i);
    % do a 5-fold cross-validation
    for j = 1:options.nCross
        modelD.cvSet = cPart{j};
        % initialize the model
        model = modelInit(data,modelD);
        % train the model
        model = modelTrain(model);
        % calculate test likelihood
        for k = modelD.cvSet
            vPre = modelPredictPiece(model,data.censor{k}.vCInt);
            vLikeBw(i) = vLikeBw(i)+sum(data.censor{k}.vCen.*log(vPre)-vPre);
            if i == 3 && options.figure == 1
                vSaveTest(k) = sum(data.censor{k}.vCen.*log(vPre)-vPre);
            end
            
        end
    end
end
fprintf('\n');
[~,nBest] = max(vLikeBw);
modelD.bandwidth = vBwSeq(nBest);
% train the model again
modelD.cvSet = [];
% initialize the model
model = modelInit(data,modelD);
% train the model
model = modelTrain(model);
%% read the timer here
nTime = toc(tStart);
% calculate the value at testing point
vT = linspace(data.XBeg,data.XEnd,options.testSample)';
model.pre = modelPredictPoint(model,vT);
            
%% test to output some data or plot
if options.figure == 1
    figure;title('LocalEM cross validation');
    plot(vBwSeq,vLikeBw,'b*-');hold on;
    if isfield(data,'fBasic') && options.dataWeight ~= 1 
        nLikeTrain = modelLikeIdeal(data,1);    % Train
        plot(vBwSeq,nLikeTrain*ones(size(vBwSeq)),'r*-'); % plot the ideal likelihood

    end
end
end

function f = calInt(model,nBeg,nEnd,options)
% calculate the integral by Simpson's rule
switch options.intType
    case 1
        np = options.intBinNum+1;
    case 2
        np = options.entBinNum+1;
end
vx = linspace(nBeg,nEnd,np)';
nh = vx(2)-vx(1);
vf = modelPredictPoint(model,vx);
f = nh/3*(vf(1)+vf(np)+4*sum(vf(2:2:np-1))...
            +2*sum(vf(3:2:np-2)));
end

function [nLike,vLike] = modelLikeIdeal(data,bflag)
% calculate the ideal likelihood
if bflag == 1 % train
    vLike = zeros(data.nTrain,1);
    for u = 1:data.nTrain
        vJ = unique([data.censor{u}.vCInt(:,1);data.censor{u}.vCInt(:,2)]);     % union
        vJl = vJ(1:end-1); % column vector
        vJr = vJ(2:end);   % column vector
        nJ = length(vJl);
        vWid = (vJr-vJl)';
        gq = GaussianQuadratic(vJl,vJr);
        vlattice = reshape(gq.lattice,gq.nQuad,nJ)';
        vInt = diag(0.5*vWid)*data.censor{u}.fCurrent(vlattice)*gq.weight';
        
        vLike(u) = sum(data.censor{u}.vCen.*log(vInt)-vInt);
    end
    nLike = sum(vLike);
else          % test
    vLike = zeros(data.nTest,1);
    for u = 1:data.nTest
        vJ = unique([data.centest{u}.vCInt(:,1);data.centest{u}.vCInt(:,2)]);     % union
        vJl = vJ(1:end-1); % column vector
        vJr = vJ(2:end);   % column vector
        nJ = length(vJl);
        vWid = (vJr-vJl)';
        gq = GaussianQuadratic(vJl,vJr);
        vlattice = reshape(gq.lattice,gq.nQuad,nJ)';
        vInt = diag(0.5*vWid)*data.centest{u}.fCurrent(vlattice)*gq.weight';
        
        vLike(u) = sum(data.centest{u}.vCen.*log(vInt)-vInt);
    end
    nLike = sum(vLike);
end
end

function vPre = modelPredictPiece(model,mCen)
% predict for given period
% mCen -- censored window
vJ = unique([mCen(:,1);mCen(:,2)]);     % union
vJl = vJ(1:end-1); % column vector
vJr = vJ(2:end);   % column vector
nJ = length(vJl);
vWid = (vJr-vJl)';
gq = GaussianQuadratic(vJl,vJr);
vlattice = gq.lattice;

[~,mK_x] = modelPredictPoint(model,vlattice);
if length(vWid) == 1
    mK = mK_x*(0.5*vWid*gq.weight');
else
    mK = mK_x*kron(diag(0.5*vWid),gq.weight');
end
vPre = model.tmpEM*mK;
if ~iscolumn(vPre)
    vPre = vPre';
end
end

function [vPre,mK_x] = modelPredictPoint(model,vT)
% predict at given points
if iscolumn(vT)
    vT = vT';
end
mLattice = repmat(vT,length(model.Jr),1);
mA = repmat(model.Jr,1,length(vT))-mLattice;
mB = repmat(model.Jl,1,length(vT))-mLattice;
mTmpK = normcdf(mA/model.bandwidth)-normcdf(mB/model.bandwidth); % J * lattice
vWid = (model.Jr-model.Jl)';

if length(vWid) > 1
    mK_x = diag(model.offsets./vWid)*mTmpK*diag(1./(model.offsets*mTmpK));
else
    error('stop');
end
vPre = model.tmpEM*mK_x;
end

function gq = GaussianQuadratic(Jl,Jr)
% Performs Gaussian Quadratic Rult
% input: J position of right-part of the intervals (Full)= [
% into row vector
if iscolumn(Jr)
    Jr = Jr';
end
if iscolumn(Jl)
    Jl = Jl';
end
pos = [-0.973906528517172,-0.865063366688985,-0.679409568299024,-0.433395394129247,...
    -0.148874338981631,0.148874338981631,0.433395394129247,0.679409568299024,...
    0.865063366688985,0.973906528517172];
vMid = 0.5*(Jr+Jl);
vWid = 0.5*(Jr-Jl);

gq.nQuad = length(pos);
gq.lattice = kron(vMid,ones(1,gq.nQuad))+kron(vWid,pos);
gq.weight = [0.066671344308688,0.149451349150581,0.219086362515982,0.269266719309996,...
    0.295524224714753,0.295524224714753,0.269266719309996,0.219086362515982,...
    0.149451349150581,0.066671344308688];
end

function model = modelInit(data,modelD)
% initialize the computation
% get information from the data
mPro = [];
vTrain = setdiff(1:data.nU,modelD.cvSet);
vJEnd = zeros(length(vTrain),1);
for i = 1:length(vTrain)
    u = vTrain(i);
    nm = max(data.censor{u}.vCInt(:,2));
    vTmp = [u*ones(data.censor{u}.nCInt,1),data.censor{u}.vCInt,data.censor{u}.vCen,nm*ones(data.censor{u}.nCInt,1)];
    mPro = [mPro;vTmp];
    vJEnd(i) = nm;
end

% precomputation
model.data = mPro;
vJLeft = model.data(:,2);
vJRight = model.data(:,3);
vJ = unique([vJLeft;vJRight]);     % union
vJl = vJ(1:end-1); % column vector
vJr = vJ(2:end);   % column vector
nJ = length(vJl);
vWid = (vJr-vJl)';
gq = GaussianQuadratic(vJl,vJr);

vlattice = gq.lattice;
vJr = vJr';
mTmp = repmat(vJEnd,1,length(vJr))- repmat(vJr,length(vJEnd),1);
voffsets = sum(mTmp >= 0,1);

% for interval with at least one event
vFlag = model.data(:,4) > 0;
vJLeft = vJLeft(vFlag);
vJRight = vJRight(vFlag);
mTmpA = repmat(vJLeft,1,length(vJr))- repmat(vJr,length(vJLeft),1);
mTmpB = repmat(vJRight,1,length(vJr))- repmat(vJr,length(vJRight),1);
mincident = (mTmpB >= 0)-(mTmpA >= 0);

vJr = vJr';
nbw = modelD.bandwidth;
mLattice = repmat(vlattice,length(vJr),1);
mA = repmat(vJr,1,length(vlattice))-mLattice;
mB = repmat(vJl,1,length(vlattice))-mLattice;
mTmpK = normcdf(mA/nbw)-normcdf(mB/nbw); % J * lattice

if length(vWid) > 1
    mK_x = diag(voffsets./vWid)*mTmpK*diag(1./(voffsets*mTmpK));
end
mK = mK_x*kron(diag(0.5*vWid),gq.weight');

model.Jr = vJr;
model.Jl = vJl;
model.incident = mincident;
model.offsets = voffsets;
model.lattice = vlattice;
model.K_x = mK_x;
model.K = mK;
model.bandwidth = nbw;
end

function model = modelTrain(model)
% use EM algorithm to train the model
vTmpM = model.data(:,4);
vTmpM = vTmpM(vTmpM>0);
vLambda = vTmpM'*model.incident*diag(1./model.offsets);
while true
    vLambdaPre = vLambda;
    mInv = diag(1./(vLambda*model.incident'));
    model.tmpEM = vTmpM'*mInv*model.incident*diag(vLambda./model.offsets);
    vLambda = model.tmpEM*model.K;
    nTmp = sum((vLambda-vLambdaPre).^2);
    if nTmp <= 1e-8
        model.lambda = model.tmpEM*model.K_x;
        break
    end
end
end



