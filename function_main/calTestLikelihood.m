function [nLike,vLike,nLikeTrain,vLikeTrain] = calTestLikelihood(data,f,vWTest,vWTrain)
% compute the likelihood for a data set
% Input:
%       data   -- should contain basic test part, see generateIntCenBook.m
%       f      -- function or struct (linear interpolation will be used if
%                 value is not found)
%                 struct:
%                   f.x: points, should be in ascending order on a grid
%                   f.y: function value
% Output:
%       nLike  -- Test likelihood over the whole set
%       vLike  -- Test likelihood for each subject
bfulllikelihood = 1; % whether to calculate the full likelihood or not

%% Input data check
if isempty(vWTest)
    vWTrain = ones(data.nTrain,1);
    vWTest = ones(data.nTest,1);
end

if ~isfield(data,'nTest') || ~isfield(data,'centest')
    error('Illegal Data Structure Input');
end
%% Compute the test likelihood by each subject
vLike = zeros(data.nTest,1);
if isempty(f)
    nEmpty = 1;
else
    nEmpty = 0;
end
for u = 1:data.nTest
    if bfulllikelihood == 1
        vLike(u) = vLike(u) - data.centest{u}.nLike;
    end
    if nEmpty == 1
        f = data.centest{u}.fCurrent;
    end
    vLike(u) = vLike(u) + sum(data.centest{u}.vCen)*log(vWTest(u));
    for ic = 1:data.centest{u}.nCInt
        if data.centest{u}.vCen(ic) > 0
            nBeg = data.centest{u}.vCInt(ic,1);
            nEnd = data.centest{u}.vCInt(ic,2);
            nInt = SimpsonRule(nBeg,nEnd,[],f);
            vLike(u) = vLike(u) + data.centest{u}.vCen(ic)*log(nInt);
        end
    end
    nBeg = min(data.centest{u}.vCInt(:,1));
    nEnd = max(data.centest{u}.vCInt(:,2));
    nInt = SimpsonRule(nBeg,nEnd,[],f);
    vLike(u) = vLike(u) - vWTest(u)*nInt;
end
nLike = sum(vLike);

%% Compute the train likelihood by each subject
vLikeTrain = zeros(data.nTrain,1);
for u = 1:data.nTrain
    if bfulllikelihood == 1
        vLikeTrain(u) = vLikeTrain(u) - data.censor{u}.nLike;
    end
    if nEmpty == 1
        f = data.censor{u}.fCurrent;
    end
    vLikeTrain(u) = vLikeTrain(u) + sum(data.censor{u}.vCen)*log(vWTrain(u));
    for ic = 1:data.censor{u}.nCInt
        if data.censor{u}.vCen(ic) > 0
            nBeg = data.censor{u}.vCInt(ic,1);
            nEnd = data.censor{u}.vCInt(ic,2);
            nInt = SimpsonRule(nBeg,nEnd,[],f);
            vLikeTrain(u) = vLikeTrain(u) + data.censor{u}.vCen(ic)*log(nInt);
        end
    end
    nBeg = min(data.censor{u}.vCInt(:,1));
    nEnd = max(data.censor{u}.vCInt(:,2));
    nInt = SimpsonRule(nBeg,nEnd,[],f);
    vLikeTrain(u) = vLikeTrain(u) - vWTrain(u)*nInt;
end
nLikeTrain = sum(vLikeTrain);
end