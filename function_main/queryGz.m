function [mGG,mDGG] = queryGz(z,gz,dgz)
%QUERYGZ Summary of this function goes here
% linear interpolation
debug = 0;
if debug
    mGG = -z;
    mDGG = ones(size(z));
else
    z = log10(z);
    mPos = 49999*(max(-11,min(11,z))+11)/22+1;
    mPos_l = floor(mPos);
    mPos_u = ceil(mPos);
    mRatio_l = (mPos_u-mPos)./(mPos_u-mPos_l);
    mRatio_u = (mPos-mPos_l)./(mPos_u-mPos_l);
    mBequal = mPos_l==mPos_u;
    mGGequal = gz(mPos_l);
    mGG = mRatio_l.*gz(mPos_l)+mRatio_u.*gz(mPos_u);
    mGG(mBequal) = mGGequal(mBequal);
    mGG(z<-11) = 0;
    if nargout == 2
        mDGGequal = dgz(mPos_l);
        mDGG = mRatio_l.*dgz(mPos_l)+mRatio_u.*dgz(mPos_u);
        mDGG(mBequal) = mDGGequal(mBequal);
        mDGG(z<-11) = 2;
    end
    
end
end

%% old version
% if z>-700
%     g= gz(floor(-z/0.001)+1);
%     dg = dgz(floor(-z/0.001)+1);
% else 
%     g = gz(700/0.001+1);
%     dg = dgz(700/0.001);