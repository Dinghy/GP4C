function plotIntCen(fig,data,options,res)
%% result from intensity censored data
figure(fig);
nLine = 1.5;
nFont = 13;
set(gca,'FontSize',nFont);
cColor = {[0,1,0],[1,0,0],[1,0,1],[0,0,1],[0,1,1]};
if nargin == 4
    model = [];
end
%% Panel Count Data 
if ~iscell(res) % Synthetic data sets with multiple experiments
    subplot('Position',[0.1,0.68,0.8,0.25]);hold on;
    vT = linspace(data.XBeg,data.XEnd,options.testInt)';
    vT = vT(1:options.testInt-1);
    % true intensity function
    vY = data.fBasic(vT);
    h1 = plot(vT,vY,'k','LineWidth',nLine);grid on;box on;
    ylim([0 max(vY)*1.5]);xlim([data.XBeg,data.XEnd]);
    % benchmark
    if ~isempty(res.model{4})
        h2 = plot(vT,res.model{4},'b','LineWidth',nLine);
    end
    % complete data
    [vMean,vHigh,vLow] = confModel(res.model{1},vT);
    h3 = plot(vT,vMean,'r','LineWidth',nLine);
    plot(vT,vHigh,'r--','LineWidth',nLine);
    plot(vT,vLow,'r--','LineWidth',nLine);
    % censored data with GP
    [vMean,vHigh,vLow] = confModel(res.model{2},vT);
    h4 = plot(vT,vMean,'m','LineWidth',nLine);
    plot(vT,vLow,'m--','LineWidth',nLine);
    plot(vT,vHigh,'m--','LineWidth',nLine);
    leg = legend([h1,h2,h3,h4],'True Intensity','Local Constant','Complete (GP)','Censored (GP)','Location','best','Orientation','horizontal');
    set(leg,'Position',[0.1,0.95,0.8,0.05]);
    ylabel('Intensity');xlabel('Time t');
    %% plot computational time
    subplot('Position',[0.1,0.45,0.8,0.15]);hold on;grid on;box on;
    nplot = 4;
    strleg = {'Complete (GP)','Benchmark','Censored (GP,L2)','Censored (GP,L3)'};
    boxplot([res.mTime(1,:)',res.mTime(4,:)',res.mTime(2,:)',res.mTime(3,:)'],...
        strleg(1:nplot));
    ylabel('Computational Time');
    %% plot MSE error
    subplot('Position',[0.1,0.25,0.8,0.15]);hold on;grid on;box on;
    boxplot([res.mError(1,:)',res.mError(4,:)',res.mErrorS(2,:)',res.mErrorS(3,:)'],strleg(1:nplot));
    ylabel('MSE');
%     switch options.dataType
%         case 1
%             ylim([0 100]);
%         case 3
%             ylim([0 50]);
%     end
    %% plot test likelihood
    subplot('Position',[0.1,0.05,0.8,0.15]);hold on;grid on;box on;
    boxplot([res.mTestLike(1,:)',res.mTestLike(4,:)',res.mTestLike(2,:)',res.mTestLike(3,:)'],strleg(1:nplot));
    ylabel('Test likelihood');
else
    for ires = 1:length(res)
        restmp = res{ires};
        ia = 0.06+(ires-1)*0.44;ib = 0.39;
        subplot('Position',[ia,0.68,ib,0.25]);hold on;
        vT = linspace(data.XBeg,data.XEnd,options.testInt)';
        vT = vT(1:options.testInt-1);
        % benchmark
        h1 = plot(vT,restmp.model{3},'b','LineWidth',nLine);
        % censored data with GP
        [vMean,vHigh,vLow] = confModel(restmp.model{2},vT);
        h2 = plot(vT,vMean,'m','LineWidth',nLine);
        plot(vT,vLow,'m--','LineWidth',nLine);
        plot(vT,vHigh,'m--','LineWidth',nLine);
        if ires == 1
            leg = legend([h1,h2],'Benchmark','Censored (GP,L2)','Location','best','Orientation','horizontal');
            set(leg,'Position',[0.1,0.95,0.8,0.05]);
        end
        ylabel('Intensity');xlabel('Time t');xlim([data.XBeg,data.XEnd]);
        %% plot computational time
        subplot('Position',[ia,0.35,ib,0.25]);hold on;grid on;box on;
        nplot = 3;
        strleg = {'Benchmark','Censored (GP,L2)','Censored (GP,L3)'};
        boxplot([restmp.mTime(3,:)',restmp.mTime(1,:)',restmp.mTime(2,:)'],...
            strleg(1:nplot));
        ylabel('Computational Time');
        %% plot test likelihood
        subplot('Position',[ia,0.05,ib,0.25]);hold on;grid on;box on;
        boxplot([restmp.mTestLike(3,:)',restmp.mTestLike(1,:)',restmp.mTestLike(2,:)'],strleg(1:nplot));
        ylabel('Test likelihood');
    end
end
%         % plot time
%         subplot('Position',[0.1,0.67,0.8,0.3]);hold on;grid on;box on;
%         vX = options.censorNum;
%         [a,b] = size(res.mFullTime);
%         vData = reshape(res.mFullTime,a*b,1);
%         [nMean,nLo,nHi] = calQuantile(vData);
%         rectangle('Position',[data.XBeg-1,nMean-nLo,vX(end)+2-data.XBeg,nHi+nLo],...
%             'FaceColor',([1,1,1]+cColor{1})/2,'EdgeColor',([1,1,1]+cColor{1})/2);
%         plot(vX,nMean*ones(size(vX)),'Color',cColor{1},'LineWidth',nLine);
%         mRes = zeros(data.nU,3);
%         for i = 1:res.nInferMax
%             for j = 1:data.nU
%                 [nm,nl,nh] = calQuantile(res.mCenTime{i}(j,:));
%                 mRes(j,:) = [nm,nl,nh]';
%             end
%             errorbar(vX,mRes(:,1),mRes(:,2),mRes(:,3),'Color',cColor{1+i},'LineWidth',nLine);
%         end
%         for i=1
%             for j = 1:data.nU
%                 [nm,nl,nh] = calQuantile(res.mBenchTime{i}(j,:));
%                 mRes(j,:) = [nm,nl,nh]';
%             end
%             errorbar(vX,mRes(:,1),mRes(:,2),mRes(:,3),'Color',cColor{4+i},'LineWidth',nLine);
%         end
%         set(gca,'XTick',vX);
%         set(gca,'FontSize',nFont);
%         xlim([vX(1)-1,vX(end)+1]);ylabel('Time(seconds)');
%         xlabel('The number of censored periods');
%         leg = legend('Complete','L1','L2','L3','B1','Location','best','Orientation','horizontal');
%         set(leg,'Position',[0.1,0.01,0.8,0.05]);
%         % plot error    
%         subplot('Position',[0.1,0.2,0.8,0.3]);hold on;grid on;box on;
%     end
    
%     if length(options.censorNum) == 1
%         nplot = 4;
%         strleg = {'Complete','L2','L3','B1'};
%         vX = options.censorNum;
%         [a,b] = size(res.mFullError);
%         mRes = zeros(options.sampleRepeat,nplot);
%         mRes(:,1) = reshape(res.mFullError,a*b,1);
%         vResCen = zeros(options.sampleRepeat,res.nInferMax-1);
%         mRes(:,2) = res.mCenError{2}(1,:)';
%         mRes(:,3) = res.mCenError{3}(1,:)';
%         if isfield(res,'mBenchError')
%             mRes(:,4) = res.mBenchError{1}(1,:)';
%         end
%         boxplot(mRes,strleg(1:nplot));
%     else
%         vX = options.censorNum;
%         [a,b] = size(res.mFullError);
%         vData = reshape(res.mFullError,a*b,1);
%         [nMean,nLo,nHi] = calQuantile(vData);
%         rectangle('Position',[data.XBeg-1,nMean-nLo,vX(end)+2-data.XBeg,nHi+nLo],...
%             'FaceColor',([1,1,1]+cColor{1})/2,'EdgeColor',([1,1,1]+cColor{1})/2);
%         plot(vX,nMean*ones(size(vX)),'Color',cColor{1},'LineWidth',nLine);
%         mRes = zeros(data.nU,3);
%         for i = 1:res.nInferMax
%             for j = 1:data.nU
%                 [nm,nl,nh] = calQuantile(res.mCenError{i}(j,:));
%                 mRes(j,:) = [nm,nl,nh]';
%             end
%             errorbar(vX,mRes(:,1),mRes(:,2),mRes(:,3),'Color',cColor{1+i},'LineWidth',nLine);
%         end
%         for i=1
%             for j = 1:data.nU
%                 [nm,nl,nh] = calQuantile(res.mBenchError{i}(j,:));
%                 mRes(j,:) = [nm,nl,nh]';
%             end
%             errorbar(vX,mRes(:,1),mRes(:,2),mRes(:,3),'Color',cColor{4+i},'LineWidth',nLine);
%         end
%         set(gca,'XTick',vX);
%         set(gca,'FontSize',nFont);
%         xlim([vX(1)-1,vX(end)+1]);
%         xlabel('The number of censored periods');
%     end
%     ylabel('mse');
% else
%     %% censoring one
%     % plot basic information
%     vT = linspace(data.XBeg,data.XEnd,200);
%     subplot('Position',[0.1,0.67,0.8,0.3]);
%     plot(vT,data.fBasic(vT),'LineWidth',nLine);grid on;box on;
%     set(gca,'FontSize',nFont);
%     ylabel('Intensity');
%     xlabel('Time t');ylim([0 ceil(max(data.fBasic(vT)))]);
%     % show performance
%     subplot('Position',[0.1,0.22,0.8,0.29]);hold on;grid on;box on;
%     vX = cellfun(@(x)x.vCInt(1),data.censor);
%     [a,b] = size(res.mFullError);
%     vData = reshape(res.mFullError,a*b,1);
%     [nMean,nLo,nHi] = calQuantile(vData);
%     
%     rectangle('Position',[data.XBeg-1,nMean-nLo,vX(end)+2-data.XBeg,nHi+nLo],...
%         'FaceColor',([1,1,1]+cColor{1})/2,'EdgeColor',([1,1,1]+cColor{1})/2);
%     plot(vX,nMean*ones(size(vX)),'Color',cColor{1},'LineWidth',nLine);
%     mRes = zeros(data.nU,3);
%     for i = 1:res.nInferMax
%         for j = 1:data.nU
%             [nm,nl,nh] = calQuantile(res.mCenError{i}(j,:));
%             mRes(j,:) = [nm,nl,nh]';
%         end
%         errorbar(vX,mRes(:,1),mRes(:,2),mRes(:,3),'Color',cColor{1+i},'LineWidth',nLine);
%     end
%     set(gca,'XTick',vX);
%     set(gca,'FontSize',nFont);
%     xlim([data.XBeg,data.XEnd]);ylabel('mse');
%     xlabel('starting time Ta of the censored period');
%     leg = legend('Complete','L1','L2','L3','Location','best','Orientation','horizontal');
%     set(leg,'Position',[0.1,0.01,0.8,0.05]);
% end
end

function vEdge = getAllCensorEdge(data)
%% get all edges from censored data
% Usage:
%         vEdge = getAllCensorEdge(data);
%         vEdge = sortrows(vEdge,[1,2]);
%         vY = linspace(-1,1,size(vEdge,1));
%         for ibar = 1:size(vEdge,1)
%             line(vEdge(ibar,:),[vY(ibar),vY(ibar)],'Color','k','LineStyle','-');
%         end
%         xlim([data.XBeg,data.XEnd]);
vEdge = [];
for u = 1:data.nU
    vEdge = [vEdge;data.censor{u}.vCInt];
end
end



function [nm,nl,nh] = calQuantile(a)
nm = mean(a);
nl = nm - quantile(a,0.25);
nh = quantile(a,0.75) - nm;
end