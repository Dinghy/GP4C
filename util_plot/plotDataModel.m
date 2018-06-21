function [fig1,fig2,fig3] = plotDataModel(data,resA,resB,options)
% Plot the training data set and testing data set and compare the test
% likelihhod of two methods
% Input:      data  --  Data set for training and testing
%             resA  --  Comparing method A
%             resB  --  Comparing method B
%          options  --  options on all parts
% Output:     fig1  --  (Train) figure handle to show the details of each subjects
%             fig2  --  (Test)  figure handle to show the details of each subjects
%             fig3  --  figure handle to show the estimated intensity

if ~isfield(options,'plotDemo');
    options.plotDemo = 0;
end

if options.plotDemo == 0
    fig1 = figure;set(fig1,'Position',[0,0,950,600]);
    vCount = [];
    for u = 1:data.nTrain
        vCount = [vCount;data.censor{u}.vCen];
    end
    for u = 1:data.nTest
        vCount = [vCount;data.centest{u}.vCen];
    end
    vCount = unique(vCount);
    nNum = length(vCount);
    cColor = mat2cell(jet(nNum),ones(size(vCount)),3);
    mapObj = containers.Map(vCount,cColor);
    nY = max(data.nTrain,data.nTest)+5;
    subplot('Position',[0.05,0.1,0.15,0.85]);hold on;box on;ylim([0 nY]);xlim([data.XBeg,data.XEnd+10]);
    for u = 1:data.nTrain
        for ic = 1:data.censor{u}.nCInt
            nBeg = data.censor{u}.vCInt(ic,1);
            nEnd = data.censor{u}.vCInt(ic,2);
            ncount = data.censor{u}.vCen(ic);
            rectangle('Position',[nBeg,u-0.25,nEnd-nBeg,0.5],'FaceColor',mapObj(ncount));
        end
    end
    title('Training Set');xlabel('Time');ylabel('Subject ID');
    subplot('Position',[0.22,0.1,0.05,0.85]);hold on;set(gca,'Visible','off'); set(gca,'XTick',[]);set(gca,'YTick',[]);
    cCount = keys(mapObj);
    nYs = nY/length(mapObj);
    for i = 1:length(mapObj)
        rectangle('Position',[0,nYs*(i-0.5)-nYs/4,1,nYs/2],'FaceColor',mapObj(cCount{i}(1)));
        text(1.5,nYs*(i-0.5),num2str(cCount{i}(1)));
    end
    ylim([0 nY]);xlim([0,2]);
    
    subplot('Position',[0.33,0.1,0.15,0.85]);hold on;box on;grid on;ylim([0 nY]);
    plot(resA.vSampTrainLike,1:data.nTrain,'r*-');
    plot(resB.vTrainLike,1:data.nTrain,'b*-');
    disp([sum(resA.vSampTrainLike),sum(resB.vTrainLike)]);
    legend('GP4C','localEM','Location','best');xlabel('Train Likelihood');ylabel('Subject ID');
    
    subplot('Position',[0.53,0.1,0.15,0.85]);hold on;box on;ylim([0 nY]);xlim([data.XBeg,data.XEnd+10]);
    for u = 1:data.nTest
        for ic = 1:data.centest{u}.nCInt
            nBeg = data.centest{u}.vCInt(ic,1);
            nEnd = data.centest{u}.vCInt(ic,2);
            ncount = data.centest{u}.vCen(ic);
            rectangle('Position',[nBeg,u-0.25,nEnd-nBeg,0.5],'FaceColor',mapObj(ncount));
        end
    end
    title('Testing Set');xlabel('Time');ylabel('Subject ID');
    subplot('Position',[0.7,0.1,0.05,0.85]);hold on;set(gca,'Visible','off'); set(gca,'XTick',[]);set(gca,'YTick',[]);
    cCount = keys(mapObj);
    for i = 1:length(mapObj)
        rectangle('Position',[0,nYs*(i-0.5)-nYs/4,1,nYs/2],'FaceColor',mapObj(cCount{i}(1)));
        text(1.5,nYs*(i-0.5),num2str(cCount{i}(1)));
    end
    ylim([0 nY]);xlim([0,2]);
    
    subplot('Position',[0.8,0.1,0.15,0.85]);hold on;box on;grid on;ylim([0 nY]);
    plot(resA.vSampLike,1:data.nTest,'r*-');
    plot(resB.vTestLike,1:data.nTest,'b*-');
    disp([sum(resA.vSampLike),sum(resB.vTestLike)]);
    legend('GP4C','localEM','Location','best');xlabel('Test Likelihood');ylabel('Subject ID');
    
    %% another plot about the intensity estimated
    fig2 = figure;set(fig2,'Position',[0,0,600,600]);
    ilinewidth = 1.5;
    vT = linspace(data.XBeg,data.XEnd,options.testSample)';
    model = resA.model;
    Kmm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.Xm);
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
    
    KnmInvKmm = Knm/Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);
    
    vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
    vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
    vMean = vMean.^2+vVar;
    
    h1 = plot(vT,vMean,'r','LineWidth',ilinewidth);hold on;
    
    plot(vT,vHigh,'r--','LineWidth',ilinewidth);
    plot(vT,vLow,'r--','LineWidth',ilinewidth);
    model = resB.model;
    h2 = plot(vT,model.pre,'m-','LineWidth',ilinewidth);hold off;
    xlim([data.XBeg,data.XEnd]);
    legend([h1,h2],'GP4C','LocalEM');
    xlabel('Time');ylabel('Intensity');
    
else % less plots
    fig1 = figure;set(fig1,'Position',[0,0,600,400]);
    vCount = [];
    for u = 1:data.nTrain
        vCount = [vCount;data.censor{u}.vCen];
    end
    for u = 1:data.nTest
        vCount = [vCount;data.centest{u}.vCen];
    end
    vCount = unique(vCount);
    nNum = length(vCount);
    mColor = jet(nNum);
    mColor(1,:) = (mColor(1,:)+[1,1,1])/2;
    cColor = mat2cell(mColor,ones(size(vCount)),3);
    mapObj = containers.Map(vCount,cColor);
    
    
    subplot('Position',[0.10,0.13,0.78,0.80]);hold on;box on;grid on;
    % subplot('Position',[0.07,0.11,0.37,0.84]);hold on;box on;
    nY = max(data.nTrain,data.nTest)+1;
    ylim([0 nY]);xlim([data.XBeg,data.XEnd+2]);
    for u = 1:data.nTrain
        for ic = 1:data.censor{u}.nCInt
            nBeg = data.censor{u}.vCInt(ic,1);
            nEnd = data.censor{u}.vCInt(ic,2);
            ncount = data.censor{u}.vCen(ic);
            rectangle('Position',[nBeg,u-0.25,nEnd-nBeg,0.5],'FaceColor',mapObj(ncount));
        end
    end
    title('Training Set in Bladder Cancer Set');xlabel('Time (Week)');ylabel('Subject (Patient) ID k');
    set(gca,'YTick',unique([0,10:5:25,6]));
    annotation('textarrow',[22 19]/(data.XEnd+2),[5.7 5.7]/nY,'String','Patient No.4','FontSize',13);
    set(gca,'FontSize',13);
    
    subplot('Position',[0.91,0.1,0.07,0.85]);hold on;set(gca,'Visible','off'); set(gca,'XTick',[]);set(gca,'YTick',[]);
    cCount = keys(mapObj);
    nYs = nY/length(mapObj);
    for i = 1:length(mapObj)
        rectangle('Position',[0,nYs*(i-0.5)-nYs/4,1,nYs/2],'FaceColor',mapObj(cCount{i}(1)));
        text(1.5,nYs*(i-0.5),num2str(cCount{i}(1)));
    end
    ylim([0 nY]);xlim([0,2]);
    set(gca,'FontSize',13);
    %% another plot about the testing data set
    fig2 = figure;set(fig2,'Position',[0,0,600,400]);
    nY = max(data.nTest)+1;
    subplot('Position',[0.03,0.12,0.07,0.80]);hold on;
    set(gca,'Visible','off'); set(gca,'XTick',[]);set(gca,'YTick',[]);
    nYs = nY/length(mapObj);
    for i = 1:length(mapObj)
        rectangle('Position',[0,nYs*(i-0.5)-nYs/4,1,nYs/2],'FaceColor',mapObj(cCount{i}(1)));
        text(1.5,nYs*(i-0.5),num2str(cCount{i}(1)));
    end
    ylim([0 nY]);xlim([0,2]);
    set(gca,'FontSize',13);
    
    subplot('Position',[0.2,0.13,0.37,0.80]);hold on;box on;grid on;
    ylim([0 nY]);xlim([data.XBeg,data.XEnd+2]);
    for u = 1:data.nTest
        for ic = 1:data.centest{u}.nCInt
            nBeg = data.centest{u}.vCInt(ic,1);
            nEnd = data.centest{u}.vCInt(ic,2);
            ncount = data.centest{u}.vCen(ic);
            rectangle('Position',[nBeg,u-0.25,nEnd-nBeg,0.5],'FaceColor',mapObj(ncount));
        end
    end
    title('Test Set in Bladder A');xlabel('Time (Week)');ylabel('Subject (Patient) ID k');
    set(gca,'FontSize',13);
    
    subplot('Position',[0.63,0.13,0.32,0.80]);hold on;box on;grid on;ylim([0 nY]);
    plot(resA.vSampLike,1:data.nTest,'r*-','LineWidth',1.5);
    plot(resB.vTestLike,1:data.nTest,'b*-','LineWidth',1.5);
    disp([sum(resA.vSampLike),sum(resB.vTestLike)]);
    if options.dataWeight == 1
        legend('GP4CW','localEM','Location','best');
    else
        legend('GP4C','localEM','Location','best');
    end
    xlabel('Test Likelihood');%ylabel('Subject ID k');
    set(gca,'FontSize',13);
    %% another plot about the intensity estimated
    fig3 = figure;set(fig3,'Position',[0,0,600,400]);
    subplot('Position',[0.11,0.13,0.85,0.85]);
    ilinewidth = 1.5;
    vT = linspace(data.XBeg,data.XEnd,options.testSample)';
    model = resA.model;
    Kmm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.Xm);
    Knm = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),vT,model.Xm);
    
    KnmInvKmm = Knm/Kmm;
    vMean = KnmInvKmm*model.var.Mu;
    vVar = sum(exp(2*model.GP.logtheta(end-1:end)))+sum(KnmInvKmm.*(Knm*(Kmm\model.var.Sigma-eye(model.m))),2);
    
    vLow = ncx2inv(0.25,1,vMean.^2./vVar).*vVar;
    vHigh = ncx2inv(0.75,1,vMean.^2./vVar).*vVar;
    vMean = vMean.^2+vVar;
    
    h1 = plot(vT,vMean,'r','LineWidth',ilinewidth);hold on;
    plot(vT,vHigh,'r--','LineWidth',ilinewidth);
    plot(vT,vLow,'r--','LineWidth',ilinewidth);
    model = resB.model;
    h2 = plot(vT,model.pre,'b-','LineWidth',ilinewidth);hold off;
    xlim([data.XBeg,data.XEnd]);
    if options.dataWeight == 1
        legend([h1,h2],'GP4CW','LocalEM');
    else
        legend([h1,h2],'GP4C','LocalEM');
    end
    grid on;box on;
    xlabel('Time (Week)');ylabel('Intensity');
    set(gca,'FontSize',13);
end
end