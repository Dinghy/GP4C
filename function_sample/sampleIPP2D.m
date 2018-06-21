function mData = sampleIPP2D(flambda,vXrange,vYrange,nRlim)
% SAMPLEIPP2D Use reject sampling to sample 2 dimensional data
% flambda: basic function
% vXrange: range in X
% vYrange: range in Y
% nRlim  : Candidate proposal function
if nargin == 3
    vX = linspace(vXrange(1),vXrange(2),100);
    vY = linspace(vYrange(1),vYrange(2),100);
    [mX,mY] = meshgrid(vX,vY);
    mF = flambda(mX,mY);
    nRlim = max(max(mF))+2;
end

vBase = [vXrange(1),vYrange(1)];
vRange = [vXrange(2)-vXrange(1),vYrange(2)-vYrange(1)];
% sampling
nS = prod(vRange);
nPoint = poissrnd(nS*nRlim);
mData = repmat(vRange,nPoint,1).*rand(nPoint,2)+repmat(vBase,nPoint,1);
% rejection
vPInt = flambda(mData(:,1),mData(:,2));
vKeepTag = rand(nPoint,1) <= vPInt./nRlim;
mData = mData(vKeepTag,:);
end

