function vF = linearmapping(vX,vTx,vTy)
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
