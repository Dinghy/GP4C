function printStat(data,res,res_rev)
% Print statistics about the data set and the result
% 1. print information about the data sets (if data is not empty)
% 2. print information about the inference result

%% print statistics about the data set
if ~isempty(data)
    nU = data.nTrain+data.nTest;
    mapObj = containers.Map('KeyType','int32','ValueType','int32');
    mWin = [];
    for u = 1:data.nTrain
        mWin = [mWin;data.censor{u}.vCInt];
        vWinAdd = (1+data.XEnd)*data.censor{u}.vCInt(:,1)+data.censor{u}.vCInt(:,2);
        for i = 1:length(vWinAdd)
            mapObj(vWinAdd(i)) = 1;
        end
    end
    for u = 1:data.nTest
        mWin = [mWin;data.centest{u}.vCInt];
        vWinAdd = (1+data.XEnd)*data.centest{u}.vCInt(:,1)+data.centest{u}.vCInt(:,2);
        for i = 1:length(vWinAdd)
            mapObj(vWinAdd(i)) = 1;
        end
    end
    vWin = unique([mWin(:,1);mWin(:,2)]);
    fprintf('Subjects %d\t Ends %d\t Interval %d\t\n',nU,length(vWin),mapObj.Count);
end
%% print statistics about the result
figure;hold on;grid on;box on;
vcolor = jet(size(res,1));
vX = 1:size(res,2);
for i = 1:size(res,1)
    if isempty(res{i,1})
        continue
    end
    %% test likelihood
    if nargin == 3 %% for a special case in realworld tests
        vTestLike = cellfun(@(x)x.nTestLike,res(i,:))+cellfun(@(x)x.nTestLike,res_rev(i,:));
        nTestLike = mean(vTestLike);
        nSDTestLike = std(vTestLike);   
    else
        vTestLike = cellfun(@(x)x.nTestLike,res(i,:));
        nTestLike = mean(vTestLike);
        nSDTestLike = std(vTestLike);
    end
    %% MISE
    if isfield(res{i,1},'nMISE')
        vMISE = cellfun(@(x)x.nMISE,res(i,:));
        nMISE = median(vMISE);
        nSDMISE = std(vMISE);
    else
        nMISE = -1;
        nSDMISE = -1;
    end
    
    %% time
    if nargin == 3
        vTime = cellfun(@(x)x.nTime,res(i,:))+cellfun(@(x)x.nTime,res_rev(i,:));
        nTime = mean(vTime);
        nSDTime = std(vTime);
    else
        vTime = cellfun(@(x)x.nTime,res(i,:));
        nTime = mean(vTime);
        nSDTime = std(vTime);
    end
    
    %% likelihood by sampling
    if isfield(res{i,1},'nSampLike')
        if nargin == 3
            % vSampLike = cellfun(@(x)x.nSampLike,res(i,:))+cellfun(@(x)x.nSampLike,res_rev(i,:));
            vSampLike = cellfun(@(x)sum(x.vSampLike),res(i,:))+cellfun(@(x)sum(x.vSampLike),res_rev(i,:));
            nSampLike = mean(vSampLike);
            nSDSampLike = std(vSampLike);
        else
            % vSampLike = cellfun(@(x)x.nSampLike,res(i,:));
            vSampLike = cellfun(@(x)sum(x.vSampLike),res(i,:));
            nSampLike = mean(vSampLike);
            nSDSampLike = std(vSampLike);
        end
        
        plot(vX',vSampLike,'Color',vcolor(i,:));
        fprintf('&MISE %.1f+%.1f\t TestLike %.1f+%.1f\t Time %.0f+%.0f\\\\\n',...
        nMISE,nSDMISE,nSampLike,nSDSampLike,nTime,nSDTime);
    else
        plot(vX',vTestLike,'Color',vcolor(i,:));
        fprintf('&MISE %.1f+%.1f\t TestLike %.1f+%.1f\t Time %.0f+%.0f\\\\\n',...
        nMISE,nSDMISE,nTestLike,nSDTestLike,nTime,nSDTime);
    end
    
    
%     if isfield(res{i,1},'kernel')
%         vKernel = cellfun(@(x)x.kernel,res(i,:));
%         nKernel = mean(vKernel);
%         nSTKernel = std(vKernel);
%         fprintf('Kernel %.4f+%.4f\n',nKernel,nSTKernel);
%     end
end
% fprintf('Test Like: GP3Int %.4f\t GP3Int(Sample) %.4f\t Bench %.4f\n',nTestLike,nSampLike,nEMLike);
% fprintf('Time: GP3Int %.4f\t Bench %.4f\n',nTime,nEMTime);

end