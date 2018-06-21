function nInt = SimpsonRule(nBeg,nEnd,nPoint,f)
% compute the integral on specific region using Simpson's Rule
% Input:
%       nBeg   -- Start point, 1 dimension
%       nEnd   -- End point, 1 dimension
%       nPoint -- Estimation point in Simpson's Rule, should be odd integer
%                 number
%       f      -- function or struct (linear interpolation will be used if
%                 value is not found)
%                 struct:
%                   f.x: points, should be in ascending order on a grid
%                   f.y: function value
% Output:
%       nInt   -- Integral over the space

%% input check
if isempty(nPoint)
    nPoint = 501;
end
if mod(nPoint,2) == 0
    nPoint = nPoint+1;
end

vX = linspace(nBeg,nEnd,nPoint)';
nh = vX(2)-vX(1);
%% compute the integral
if isa(f,'function_handle') % is a function handle
    vF = f(vX);
    nInt = nh/3*(vF(1)+vF(nPoint)+4*sum(vF(2:2:nPoint-1))+2*sum(vF(3:2:nPoint-2)));

elseif isstruct(f) && isfield(f,'x') && isfield(f,'y') % is a struct
    vTx = f.x;
    vTy = f.y;
    % map vX to vF according to finer grid vTx to vTy
    if ~iscolumn(vTy)
        vTy = vTy';
    end
    nMin = min(vTx);
    nMax = max(vTx);
    nStep = vTx(2)-vTx(1);
    vPos = (min(max(nMin,vX),nMax)-nMin)/nStep+1;
    vPos_l = floor(vPos);
    vPos_u = min(length(vTy),ceil(vPos));
    vRatio_l = (vPos_u-vPos)./(vPos_u-vPos_l);
    vRatio_u = (vPos-vPos_l)./(vPos_u-vPos_l);
    
    vEqual = vPos_l == vPos_u;
    vFequal = vTy(vPos_l);
    vF = (vRatio_l.*vTy(vPos_l)+vRatio_u.*vTy(vPos_u));
    vF(vEqual) = vFequal(vEqual);
    vF(vX<nMin) = vTy(1);
    vF(vX>nMax) = vTy(end);
    nInt = nh/3*(vF(1)+vF(nPoint)+4*sum(vF(2:2:nPoint-1))+2*sum(vF(3:2:nPoint-2)));
else
    error('Illegal f Input');
end
end