function [model,nTime,vELBO] = varTrainGP3Book(data,options)
%% GP3 Model
% determine whether to plot the details of the training

%% Check the field of the training
if ~isfield(options,'detailPlot')
    options.detailPlot = 0;
end
if ~isfield(options,'approxType')
    options.approxType = 0;
    options.approxEps = 0;
end
%% determine whether to do the plots
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
%% training
load('gzdgz.mat');
model = varInit(data,options);
model.gz=gz;
model.dgz=dgz;
nTime = 0;
[LB,UB] = optConstraint(model);
nELBOpre = -Inf;
iIter = 1;
modelD = preKernel(model);
vELBO = zeros(1,options.maxIter);
tStart = tic;
while iIter <= options.maxIter 
    if options.detailPlot == 1
        para.Mu(iIter) = sum(abs(model.var.Mu));
        para.L(iIter) = sum(sum(abs(model.var.L)));
        para.theta(:,iIter) = model.GP.logtheta;
    end
    modelD = preVar(data,model,modelD);
    if options.debug == 1
        V = extractMuSigma(model);
        diff = checkGrad('calDVar',V,data,model,modelD,options);
        disp(diff);
        V = extractTheta(model,options);
        diff = checkGrad('calDKer',V,data,model,modelD,options);
        disp(diff);
        V = model.GP.logtheta;
        diff = checkMatrixGrad('MatrixPhi',V,model);
        disp(diff);
        return;
    end
    if options.figure == 1
        plotResultF(model,modelD,options,data,fig);
    end
    %% optimize mu,sigma
    vVar=extractMuSigma(model);
    myObj = @(VV)calDVar(VV,data,model,modelD,options);
    [V,nELBO] = minConf_TMP(myObj,vVar,LB,UB,struct('verbose',0,'maxIter',10));
    % return the parameter
    model = returnMuSigma(model,V);
    modelD = preVar(data,model,modelD); % update
    %% optimize hyper
    if options.hyper == 1
        vTheta = extractTheta(model,options);
        myObj = @(VV)calDKer(VV,data,model,modelD,options);
        [V,~] = minConf_TMP(myObj,vTheta,...
            model.GP.cons(1,1:length(vTheta))',model.GP.cons(2,1:length(vTheta))',struct('verbose',0,'maxIter',10));
        model = returnTheta(model,V,options);
        modelD = preKernel(model);
    end
    %% break if finished
    vELBO(iIter) = nELBO;
    if abs(nELBO-nELBOpre)<=1e-4*abs(nELBO)
        nTime = toc(tStart);
        vELBO = vELBO(1:iIter);
        iIter = iIter + 1;
        break
    else
        nELBOpre = nELBO;
        iIter = iIter + 1;
    end
end
%% Plot the detail
if options.detailPlot == 1
    plotDetail(fig2,para,iIter-1);
end
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

function model = varInit(data,options)
%% Initialization part for the model
model.const.C = 0.5772156649;
model.U = data.nU;
model.m = options.pseudo;
model.Tstart = data.XBeg;              % time start
model.Tend = data.XEnd;                % time end
model.NumData = zeros(1,data.nU);
model.wModel = options.wModel;
for u = 1:data.nU
    model.NumData(u) = size(data.train{u}.X,1);  % number of points in each dataset
end
model.cov = {'covSum',{'covSEard','covNoise'}};

model.prior.gbase = sqrt(mean(model.NumData)/prod(model.Tend-model.Tstart));
model.prior.g = model.prior.gbase*ones(model.m,1);
D = size(data.XBeg,2);

if D == 1
    model.Xm = linspace(model.Tstart,model.Tend,model.m)';
    model.plotXm = linspace(model.Tstart,model.Tend,100)';
elseif D == 2
    nPseudo = floor(sqrt(model.m));
    model.m = nPseudo^2;
    nPadX = (data.XEnd(1)-data.XBeg(1))/10;
    nPadY = (data.XEnd(2)-data.XBeg(2))/10;
    vX = linspace(data.XBeg(1)+nPadX,data.XEnd(1)-nPadX,nPseudo);
    vY = linspace(data.XBeg(2)+nPadY,data.XEnd(2)-nPadY,nPseudo);
    [mX,mY] = meshgrid(vX,vY);
    model.Xm = [reshape(mX,model.m,1),reshape(mY,model.m,1)];
    
    model.plot.m = 20;
    vX = linspace(data.XBeg(1),data.XEnd(1),model.plot.m);
    vY = linspace(data.XBeg(2),data.XEnd(2),model.plot.m);
    [mX,mY] = meshgrid(vX,vY);
    model.plot.X = mX;
    model.plot.Y = mY;
else
    error('Illegal Input');
end

if options.hyper == 0
    model.GP.logtheta = [1.18956749071803*ones(D,1);0.00748204591068048;-6.90775527898214];
else
    model.GP.logtheta = zeros(D+2,1);
    dd = log((model.Tend-model.Tstart)')/2;
    model.GP.logtheta(1:D) = dd;
    model.GP.logtheta(D+1) = 0.5*log(0.5*model.prior.gbase^2);%0.5*log(0.5*model.prior.gbase^2);
    model.GP.logtheta(D+2) = log(1e-3);
end

model.GP.cons(1,1:D) = -inf;%log((model.Tend-model.Tstart)')*0.01;
model.GP.cons(2,1:D) = inf;%log((model.Tend-model.Tstart)')*10;
model.GP.cons(1,D+1) = -inf;%log(1e-3*model.prior.gbase);%-inf;%
model.GP.cons(2,D+1) = inf;%4*log(model.prior.gbase);%inf;%
model.GP.cons(1,D+2) = log(1e-2);
model.GP.cons(2,D+2) = log(1e-5);

% Version 2: prior for q(m,S) Sophisticated 
vPos = [];
for u = 1:data.nTrain
    vPos = [vPos;data.train{u}.X];
end
vEdge = model.Xm;nStep = model.Xm(2)-model.Xm(1);
vEdge = vEdge-nStep/2;vEdge = [vEdge;model.Xm(end)+nStep/2];
vCount = histc(vPos,vEdge)/data.nTrain/nStep;
vCount = vCount(1:end-1);
vCount(1) = vCount(1)*2;vCount(end) = vCount(end)*2;
model.var.Mu = sqrt(0.5*vCount).*(0.1+rand(model.m,1)); 

% Version 1: 
% model.var.Mu = sqrt(0.5*model.prior.gbase^2)*(0.1+rand(model.m,1)); 

model.var.Sigma = feval(model.cov{:},model.GP.logtheta,model.Xm);
model.var.L = chol(model.var.Sigma,'lower');
end

function norm = checkMatrixGrad(objFunc,vInput,model) 
% function to check whether gradient caculation is correct
% 
nEps=1e-5;
norm = 0;
[~,df]=feval(objFunc,vInput,model);

for i=1:length(vInput)     
  vInput(i) = vInput(i) - 0.5*nEps;
  E0 = feval(objFunc,vInput,model);
  vInput(i) = vInput(i) + nEps;
  E1 =feval(objFunc,vInput,model);
  vInput(i) = vInput(i) - 0.5*nEps;
  diff_grad = (E1 - E0)/nEps;
  nPlus = sum(sum(abs(diff_grad - df{i})));
  norm = norm + nPlus;
  
  fprintf('%d: Diff %f\n',i,nPlus);
end
end 

function norm = checkGrad(objFunc,vInput,varargin) 
% DEBUG: function to check whether gradient caculation is correct
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

function vV = extractTheta(model,options)
%% extract theta
switch options.hyperOpt
    case 0
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta(1:end-1);
        else
            vV = model.GP.logtheta(1:end-1)';
        end
    case 1
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta(1:end-2);
        else
            vV = model.GP.logtheta(1:end-2)';
        end
    case 2
        if iscolumn(model.GP.logtheta)
            vV = model.GP.logtheta;
        else
            vV = model.GP.logtheta';
        end
end
end

function model = returnTheta(model,V,options)
%% return theta
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
nStart =nStart+M;
mtmp(mtmp==1)=V(nStart:nStart+M2-1);
model.var.L = mtmp;
model.var.Sigma = mtmp*mtmp';
end

function modelD = preKernel(model,modelD)
%  PRE-CALCULATION FOR Parameters that are only related to kernel
%  hyper-parameters
if nargin == 1 % initialization
    modelD.numPar = length(model.GP.logtheta);
end
modelD.diagKnn = sum(exp(2*model.GP.logtheta(end-1:end)));
modelD.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xm);
modelD.KmmDx = cell(1,modelD.numPar);
for k=1:modelD.numPar
    modelD.KmmDx{k} = feval(model.cov{:}, model.GP.logtheta, model.Xm, [], k);
end
% calculate integral
[mPhi,DmPhidx] = MatrixPhi(model.GP.logtheta, model);
modelD.mPhi = mPhi;
modelD.invKmmPhi = modelD.Kmm\mPhi;
modelD.DmPhidx = DmPhidx;
end

function modelD = preVar(data,model,modelD)
%  PRE-CALCULATION FOR Parameters used in update for m,S
D =length(model.Tend);

modelD.InvKmmMu = modelD.Kmm\model.var.Mu;
modelD.InvKmmSigma = modelD.Kmm\model.var.Sigma;
modelD.InvKmmL = modelD.Kmm\model.var.L;
nGamma = exp(2*model.GP.logtheta(D+1));
InvKmmSigmaMu=modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu';
modelD.vfInt=nGamma*prod(model.Tend-model.Tstart)-trace(modelD.invKmmPhi)+...
    trace(modelD.invKmmPhi*InvKmmSigmaMu);

% iterate for per batch (Book)
modelD.Knm = cell(1,data.nU);
modelD.DknmDx = cell(1,data.nU);
modelD.KnmInvKmm = cell(1,data.nU);
for u = 1:data.nU
    if ~isempty(data.train{u}.X) 
        sRes = calCovSEard(model.GP.logtheta(1:end-1),data.train{u}.X,model.Xm);
        modelD.Knm{u} = sRes.Knm;
        modelD.DknmDx{u} = sRes.DknmDx;
        modelD.KnmInvKmm{u} = sRes.Knm/modelD.Kmm;
    end
end
end

function [f,df] = calDVar(V,data,model,modelD,options)
% calculate the derivative and the current likelihood
model = returnMuSigma(model,V);
modelD = preVar(data,model,modelD);
modelD.f = 0;
modelD.var.Mu = zeros(size(model.var.Mu));
modelD.var.L = zeros(size(model.var.L));

for u=1:data.nU
    if isempty(data.train{u}.X)
        continue;
    end
    vMean = modelD.Knm{u}*modelD.InvKmmMu;
    vVar = modelD.diagKnn+sum(modelD.KnmInvKmm{u}.*(modelD.Knm{u}*(modelD.InvKmmSigma-eye(model.m))),2);
    switch options.approxType
        case 0
            vMean2DVar = vMean.^2/2./vVar;
            [vGG,vDGG] = queryGz(vMean2DVar,model.gz,model.dgz);
            mTmp_basef = 0.5*exp(-model.const.C)*exp(-vGG).*vVar;
            modelD.f = modelD.f+sum(log(mTmp_basef));
            vMutmp = vDGG.*vMean./vVar;
            vLtmp = (-vDGG.*vMean2DVar+1)./vVar;
        case 1 % mu^2
            modelD.f = modelD.f + sum(log(vMean.^2+options.approxEps)); % to guarantee it is not negative
            vMutmp = 2*vMean./(vMean.^2+options.approxEps);
            vLtmp = zeros(size(vMutmp));
        case 2 % mu^2+0.3*sigma^2
            modelD.f = modelD.f + sum(log(vMean.^2+0.3*vVar+options.approxEps));
            vMutmp = 2*vMean./(vMean.^2+0.3*vVar+options.approxEps);
            vLtmp = 0.3./(vMean.^2+0.3*vVar+options.approxEps);
        case 3 % mu^2+sigma^2
            modelD.f = modelD.f + sum(log(vMean.^2+vVar+options.approxEps));
            vMutmp = 2*vMean./(vMean.^2+vVar+options.approxEps);
            vLtmp = 1./(vMean.^2+vVar+options.approxEps);
    end
    mLMulti=repmat(vLtmp,1,model.m)';
    modelD.var.Mu = modelD.var.Mu+modelD.KnmInvKmm{u}'*vMutmp;
    modelD.var.L = modelD.var.L+(mLMulti.*modelD.KnmInvKmm{u}')*modelD.KnmInvKmm{u};
end
modelD.var.L = 2*modelD.var.L*model.var.L;
%% add the integral part and tune the learning rate rho
modelD.f = modelD.f-data.nU*modelD.vfInt;
modelD.var.Mu = modelD.var.Mu-2*data.nU*modelD.invKmmPhi*modelD.InvKmmMu;
modelD.var.L = modelD.var.L-2*data.nU*modelD.invKmmPhi*modelD.InvKmmL;
%% KL divergence and gradient
InvKmmMuG = modelD.Kmm\(model.var.Mu-model.prior.g);
InvKmmSigmaMuG=modelD.InvKmmSigma+InvKmmMuG*(model.var.Mu-model.prior.g)';
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
modelD = preKernel(model,modelD);
modelD = preVar(data,model,modelD);
%% initializing and calculation
modelD.f=0;
modelD.GP.logtheta = zeros(size(model.GP.logtheta));
modelD.InvKmmSigmaInvKmm = modelD.InvKmmSigma/modelD.Kmm;

for u = 1:data.nU
    if isempty(data.train{u}.X)
        continue;
    end
    %% point value and derivatives
    vMean = modelD.Knm{u}*modelD.InvKmmMu;
    vVar = modelD.diagKnn+...
        sum(modelD.KnmInvKmm{u}.*(modelD.Knm{u}*(modelD.InvKmmSigma-eye(model.m))),2);
    vMean2DVar = vMean.^2/2./vVar;
    switch options.approxType
        case 0
            [vGG,vDGG] = queryGz(vMean2DVar,model.gz,model.dgz);
            mTmp_basef = 0.5*exp(-model.const.C)*exp(-vGG).*vVar;
            modelD.f = modelD.f+sum(log(mTmp_basef));
            vMutmp = vDGG./vVar;
            vLtmp = 2*(-vDGG.*vMean2DVar+1)./vVar;
        case 1 % mu^2
            modelD.f = modelD.f + sum(log(vMean.^2+options.approxEps)); % to guarantee it is not negative
            vMutmp = 2./(vMean.^2+options.approxEps);
            vLtmp = zeros(size(vMutmp));
        case 2 % mu^2+0.3*sigma^2
            modelD.f = modelD.f + sum(log(vMean.^2+0.3*vVar+options.approxEps));
            vMutmp = 2./(vMean.^2+0.3*vVar+options.approxEps);
            vLtmp = 0.3./(vMean.^2+0.3*vVar+options.approxEps)*2;
        case 3 % mu^2+sigma^2
            modelD.f = modelD.f + sum(log(vMean.^2+vVar+options.approxEps));
            vMutmp = 2./(vMean.^2+vVar+options.approxEps);
            vLtmp = 1./(vMean.^2+vVar+options.approxEps)*2;
    end
    
    DbDxtmp = 2*modelD.Knm{u}*modelD.InvKmmSigmaInvKmm-modelD.KnmInvKmm{u};
    for k=1:modelD.numPar
        if k>=modelD.numPar-1
            DknnDx=2*exp(2*model.GP.logtheta(k))*ones(model.NumData(u),1);
        else % except last parameter noise
            DknnDx=zeros(model.NumData(u),1);
        end
        DknmDx = modelD.DknmDx{u}{k};
        DknmkmmInvDx = (-modelD.KnmInvKmm{u}*modelD.KmmDx{k} + DknmDx);
        DbDx = DknnDx -sum(modelD.KnmInvKmm{u}.*DknmDx,2)+...
            sum(DknmkmmInvDx.*DbDxtmp,2);
        DaDx = DknmkmmInvDx*modelD.InvKmmMu;
        modelD.GP.logtheta(k)=modelD.GP.logtheta(k)+(vMean.*DaDx)'*vMutmp+DbDx'*vLtmp/2;
    end
end
%% Integration over Tlim
nGamma = exp(2*model.GP.logtheta(length(model.Tend)+1));
InvKmmSigmaMu=modelD.InvKmmSigma+modelD.InvKmmMu*model.var.Mu';
modelD.f=modelD.f-data.nU*modelD.vfInt;
mtmp_a = (eye(model.m)-InvKmmSigmaMu)/modelD.Kmm;
mtmp_b = (2*modelD.invKmmPhi*InvKmmSigmaMu-modelD.invKmmPhi)/modelD.Kmm;
modelD.GP.logtheta = modelD.GP.logtheta+...
    data.nU*cellfun(@(x,y)trace(mtmp_a*x+mtmp_b*y),modelD.DmPhidx,modelD.KmmDx)';
if length(W) > length(model.Tend)
    modelD.GP.logtheta(end-1) = modelD.GP.logtheta(end-1)-data.nU*nGamma*2*prod(model.Tend-model.Tstart);
end    
%% KL divergence and gradient
InvKmmMuG = modelD.Kmm\(model.var.Mu-model.prior.g);
InvKmmSigmaMuG=modelD.InvKmmSigma+InvKmmMuG*(model.var.Mu-model.prior.g)';
modelD.f = modelD.f+0.5*(model.m-trace(InvKmmSigmaMuG))...
    +sum(log(diag(model.var.L)))-0.5*logdet(modelD.Kmm);
mtmp = InvKmmSigmaMuG/modelD.Kmm;
modelD.GP.logtheta = modelD.GP.logtheta-...
    0.5*cellfun(@(x)trace(modelD.Kmm\x-mtmp*x),modelD.KmmDx)';
%% extract the variable
f = - modelD.f;
df = - extractTheta(modelD,options);
end

%% plot latent function and confidence interval
function plotResultF(model,modelD,options,data,fig)
ilinewidth = 2;
if length(data.XBeg) == 2
    figure(fig);
    subplot('Position',[0.1,0.15,0.37,0.53]);
    if isfield(data,'fBasic')
        mF = data.fBasic(model.plot.X,model.plot.Y);
        imagesc(flipud(mF));title('True Intensity');
    end
    subplot('Position',[0.53,0.15,0.37,0.53]);
    mPlotXm = [reshape(model.plot.X,model.plot.m^2,1),reshape(model.plot.Y,model.plot.m^2,1)];
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),mPlotXm,model.Xm);
    KnmInvKmm = Knm/modelD.Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    mPre = reshape(vMean.^2,model.plot.m,model.plot.m);
    imagesc(flipud(mPre));title('Inferred Intensity');
    colorbar('Position',[0.92,0.15,0.03,0.53]);
    
elseif length(data.XBeg) == 1
    vT = model.plotXm;
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
    KnmInvKmm = Knm/modelD.Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    vVar = modelD.diagKnn+sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma-eye(model.m))),2);
    
    figure(fig);
    subplot('Position',[0.1,0.2,0.8,0.68]);grid on;box on;
    vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
    vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
    vMean = vMean.^2+vVar;
    
    if options.dataWeight == 0
        plot(vT,data.fBasic(vT),'k','LineWidth',ilinewidth);hold on;
        plot(vT,vMean,'r','LineWidth',ilinewidth);
    else
        plot(vT,vMean,'r','LineWidth',ilinewidth);hold on;
    end
    plot(vT,vHigh,'r--','LineWidth',ilinewidth);
    plot(vT,vLow,'r--','LineWidth',ilinewidth);hold off;
    xlim([data.XBeg,data.XEnd]);
    
    set(gca,'FontSize',12);

end
drawnow;
end
