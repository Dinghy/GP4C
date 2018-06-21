function [vMean,vHigh,vLow] = confModel(model,vT)
% return confidence interval
if iscell(model.GP.logtheta)
    hyp = model.GP.logtheta{1};
    Mu = model.var.Mu{1};
    Sigma = model.var.Sigma{1};
else
    hyp = model.GP.logtheta;
    Mu = model.var.Mu;
    Sigma = model.var.Sigma;
end
Knm = feval(model.cov{2}{1},hyp(1:end-1),vT,model.Xm);
modelD.Kmm = feval(model.cov{:},hyp,model.Xm);
modelD.diagKnn = sum(exp(2*hyp(end-1:end)));
modelD.InvKmmSigma = modelD.Kmm\Sigma;
KnmInvKmm = Knm/modelD.Kmm;
vMean = KnmInvKmm*Mu;
vVar = modelD.diagKnn+sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma-eye(model.m))),2);
vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
vMean = vMean.^2+vVar;
end