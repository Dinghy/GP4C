%% plot latent function and confidence interval
function plotIntensity(model,data,fig)
ilinewidth = 2;
if isfield(model,'flam')
    figure(fig);
    vT = model.plotXm;
    subplot('Position',[0.1,0.2,0.8,0.34]);hold on;grid on;box on;
    plot(vT,model.flam(vT),'m','LineWidth',ilinewidth);
else
    vT = model.plotXm;
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
    modelD.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xm);
    modelD.diagKnn = sum(exp(2*model.GP.logtheta(end-1:end)));
    modelD.InvKmmSigma = modelD.Kmm\model.var.Sigma;
    KnmInvKmm = Knm/modelD.Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    vVar = modelD.diagKnn+sum(KnmInvKmm.*(Knm*(modelD.InvKmmSigma-eye(model.m))),2);
    vMean2DVar = vMean.^2/2./vVar;
    vGG = queryGz(vMean2DVar,model.gz,model.dgz);
    vSig = sqrt(vVar);
    
    figure(fig);
    
    subplot('Position',[0.1,0.62,0.8,0.34]);
    semilogy(vT,0.5*exp(-model.const.C)*exp(-vGG).*vVar,'r');hold on;
    semilogy(vT,vMean.^2,'b');
    semilogy(vT,vMean.^2+0.3*vVar,'m');
    semilogy(vT,vMean.^2+vVar,'g');hold off;grid on;
    xlim([data.XBeg,data.XEnd]);
    
    subplot('Position',[0.1,0.2,0.8,0.34]);grid on;box on;
    if isfield(data,'fMean')
        vMeanTrue = data.fMean{model.nCurU}(vT);
        vSigTrue = sqrt(data.fvar{model.nCurU}(vT));
        h1 = plot(vT,vMeanTrue,'b','LineWidth',ilinewidth);hold on;
        plot(vT,vMeanTrue-2*vSigTrue,'b--','LineWidth',ilinewidth);
        plot(vT,vMeanTrue+2*vSigTrue,'b--','LineWidth',ilinewidth);
        
        h2 = plot(vT,vMean,'r','LineWidth',ilinewidth);
        plot(vT,vMean-2*vSig,'r--','LineWidth',ilinewidth);
        plot(vT,vMean+2*vSig,'r--','LineWidth',ilinewidth);
        h3 = plot(data.vXTrain,data.vFTrain,'b*');hold off;
        xlim([data.XBeg,data.XEnd]);
        legend([h1,h2,h3],'True','Inferred','Points','Location','best');
    else
        plot(vT,vMean,'r','LineWidth',ilinewidth);%hold on;
        plot(vT,vMean-2*vSig,'r--','LineWidth',ilinewidth);
        plot(vT,vMean+2*vSig,'r--','LineWidth',ilinewidth);hold off;
        xlim([data.XBeg,data.XEnd]);
    end
end
set(gca,'FontSize',12);
% plot points the bottom plot
subplot('Position',[0.1,0.05,0.8,0.05]);
for ibar=1:length(data.train{model.nCurU}.X)
    iX=data.train{model.nCurU}.X(ibar);
    line([iX,iX],[-1,1],'Color','k','LineStyle','-');
end
set(gca,'XTick',[]);set(gca,'YTickLabel',[]);ylabel('Train');
xlim([data.XBeg,data.XEnd]);
if isfield(data,'test')
    subplot('Position',[0.1,0.09,0.8,0.05]);
    for ibar=1:length(data.test{model.nCurU}.X)
        iX=data.test{model.nCurU}.X(ibar);
        line([iX,iX],[-1,1],'Color','k','LineStyle','-');
    end
    set(gca,'XTick',[]);set(gca,'YTickLabel',[]);ylabel('Test');
    xlim([data.XBeg,data.XEnd]);
end
drawnow;

end