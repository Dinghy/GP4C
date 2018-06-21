function model = kernelSmoothing(data,options)
%% kernel smoothing method 
% A Kernel Method for Smoothing Point Process Data
% By PETER DIGGLE, 1984
if ~isfield(options,'endCorrect')
    options.endCorrect = 0;
end
    
model.nCurU = 1;
model.Tstart = data.XBeg;
model.Tend = data.XEnd;
model.plotXm = linspace(model.Tstart,model.Tend,100)';
%% learn a proper window size h Gaussian kernel with various bandwidth
model.para.h = determineH(model,data);
model.fk = @(x)1/sqrt(2*pi)/model.para.h*exp(-x.^2/model.para.h^2/2);
%% return the estimated lambda
vX = data.train{model.nCurU}.X;
nX = length(vX);
if options.endCorrect == 1
    model.flam = @(x)sum(model.fk(repmat(x,1,nX)-repmat(vX',length(x),1)),2)./...
        (normcdf((data.XEnd-x)/model.para.h)-normcdf((data.XBeg-x)/model.para.h));
else
    model.flam = @(x)sum(model.fk(repmat(x,1,nX)-repmat(vX',length(x),1)),2);
end

end

function h = determineH(model,data)
% determine the hyper-parameter
% by maximising the leave-one-out training objective
vh = linspace(max(0.05,data.XBeg),data.XEnd/4,60);
vres = zeros(size(vh));
vX = data.train{model.nCurU}.X;
for i = 1:length(vh)
    nh = vh(i);
    fk = @(x)1/sqrt(2*pi)/nh*exp(-x.^2/nh^2/2);
    mTmp = fk(repmat(vX,1,length(vX))-repmat(vX',length(vX),1));
    mTmp(1:length(vX)+1:end) = zeros(1,length(vX));
    vres(i) = sum(log(sum(mTmp,2)));
end
[~,im] = max(vres);
h = vh(im);
end