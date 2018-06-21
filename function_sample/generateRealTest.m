function [dataFull,strfileName] = generateRealTest(options)
switch options.dataType 
    case 1
        dataFull = generateNausea();
        strfileName = 'Nausea';
        ncell = 2;
    case 2
        dataFull = generateBladder();
        strfileName = 'Bladder';
        ncell = 2;
    case 3
        dataFull = generateSkin();
        strfileName = 'Skin';
        ncell = 4;
end
% calculate the constant term
for c=1:ncell
    for i = 1:dataFull{c}.nTrain
        dataFull{c}.censor{i}.nLike = computeConst(dataFull{c}.censor{i}.vCen);
    end
    
    for i = 1:dataFull{c}.nTest
        dataFull{c}.centest{i}.nLike = computeConst(dataFull{c}.centest{i}.vCen);
    end
end
end

function dataFull = generateNausea()
% small function to generate Nausea data set
t = readtable('Nausea.txt','Delimiter',' ','ReadVariableNames',false);
[nrow,ncol] = size(t);
dataDose.nU = 65;
dataDose.nTrain = 0;
dataDose.nTest = 0;
dataDose.censor = cell(1,65);      % for training
dataDose.centest = cell(1,65);     % for testing
dataDose.XBeg = 0;
dataDose.XEnd = 0;
vRandA = randperm(dataDose.nU);
vRandA = vRandA(1:round(dataDose.nU/2));

dataPlacebo.nU = nrow-dataDose.nU;
dataPlacebo.nTrain = 0;
dataPlacebo.nTest = 0;
dataPlacebo.censor = cell(1,nrow-65);
dataPlacebo.censor = cell(1,nrow-65);
dataPlacebo.XBeg = 0;
dataPlacebo.XEnd = 0;
vRandB = dataDose.nU+randperm(dataPlacebo.nU);
vRandB = vRandB(1:round(dataPlacebo.nU/2));

nEndMax = 0;
nZeroRow = 0;
for i = 1:nrow            % each patient
    u = t{i,1};           % user id
    celltmp.vCInt = [];
    celltmp.nCInt = 0;
    celltmp.vLen = [];
    celltmp.vCen = [];
    celltmp.X = [];
    
    nBeg = 0;
    for j = 2:ncol        % each row
        if mod(j,2) == 0  % interval
            if isa(t{i,j},'cell')   % is a string in format
                tmp = t{i,j}{1};
                if strcmp(tmp,'.')==1 % stop at here
                    break % for this patient only
                else
                    nEnd = str2double(tmp);
                end
            else
                nEnd = t{i,j};
            end
            nEndMax = max(nEndMax,nEnd);
        else              % interval count
            if isa(t{i,j},'cell')            % is a string
                tmp = t{i,j}{1};
                nCount = str2double(tmp); % for this patient only
            else
                nCount = t{i,j};
            end
            celltmp.vCInt = [celltmp.vCInt;nBeg,nEnd];
            celltmp.nCInt = celltmp.nCInt + 1;
            celltmp.vLen = [celltmp.vLen;nEnd-nBeg];
            celltmp.vCen = [celltmp.vCen;nCount];
            nBeg = nEnd;
        end
    end
    if sum(celltmp.vCen) == 0
        nZeroRow = nZeroRow+1;
    end
    if u > 65  % Placebo Group
        if any(vRandB == u)
            dataPlacebo.nTrain = dataPlacebo.nTrain + 1;
            dataPlacebo.censor{dataPlacebo.nTrain} = celltmp;
        else
            dataPlacebo.nTest = dataPlacebo.nTest + 1;
            dataPlacebo.centest{dataPlacebo.nTest} = celltmp;
        end 
    else       % High Dose Group
        if any(vRandA == u)
            dataDose.nTrain = dataDose.nTrain + 1;
            dataDose.censor{dataDose.nTrain} = celltmp;
        else
            dataDose.nTest = dataDose.nTest + 1;
            dataDose.centest{dataDose.nTest} = celltmp;
        end
    end
end
% fprintf('Total %d\t Zero %d\n',nrow,nZeroRow);
dataDose.XEnd = nEndMax;
dataPlacebo.XEnd = nEndMax;
% save('dataNausea.mat','dataDose','dataPlacebo');
dataFull = cell(1,2);
dataFull{1} = ShrinkSet(dataDose);
dataFull{2} = ShrinkSet(dataPlacebo);
end

function dataFull = generateBladder()
% small function to generate Bladder data set
t = readtable('Bladder.txt','Delimiter',' ','ReadVariableNames',false);
tc = readtable('BladderContinue.txt','Delimiter',' ','ReadVariableNames',false);
[nrow,ncol] = size(t);
[~,ncolt] = size(tc);

dataPlacebo.nU = 47;
dataPlacebo.nTrain = 0;dataPlacebo.nTest = 0;
dataPlacebo.censor = cell(1,47);
dataPlacebo.censor = cell(1,47);
dataPlacebo.XBeg = 0;dataPlacebo.XEnd = 53;
vRandB = randperm(dataPlacebo.nU);
vRandB = vRandB(1:round(dataPlacebo.nU/2));

dataDose.nU = nrow-47;
dataDose.nTrain = 0;dataDose.nTest = 0;
dataDose.censor = cell(1,nrow-47);   % for training
dataDose.centest = cell(1,nrow-47);     % for testing
dataDose.XBeg = 0;dataDose.XEnd = 53;
vRandA = dataPlacebo.nU+randperm(dataDose.nU);
vRandA = vRandA(1:round(dataDose.nU/2));


for i = 1:nrow
    nBeg = 0;nEnd = 1;
    u = t{i,1};           % user id
    
    celltmp.vCInt = [];
    celltmp.nCInt = 0;
    celltmp.vLen = [];
    celltmp.vCen = [];
    celltmp.X = [];
    for j = 4:ncol  % first file from the fourth column
        if isa(t{i,j},'cell') % contains string
            tmp = t{i,j}{1};
            if strcmp(tmp,'.') ~= 1 % visit here
                nCount = str2double(tmp);
                celltmp.vCInt = [celltmp.vCInt;nBeg,nEnd];
                celltmp.nCInt = celltmp.nCInt + 1;
                celltmp.vLen = [celltmp.vLen;nEnd-nBeg];
                celltmp.vCen = [celltmp.vCen;nCount];
                nBeg = nEnd;
            end
            nEnd = nEnd + 1;
        else
            nCount = t{i,j};
            celltmp.vCInt = [celltmp.vCInt;nBeg,nEnd];
            celltmp.nCInt = celltmp.nCInt + 1;
            celltmp.vLen = [celltmp.vLen;nEnd-nBeg];
            celltmp.vCen = [celltmp.vCen;nCount];
            nBeg = nEnd;
        end
    end
    for j = 2:ncolt % second file
        if isa(tc{i,j},'cell') % contains string
            tmp = tc{i,j}{1};
            if strcmp(tmp,'.')~=1 % visit here 
                nCount = str2double(tmp);
                celltmp.vCInt = [celltmp.vCInt;nBeg,nEnd];
                celltmp.nCInt = celltmp.nCInt + 1;
                celltmp.vLen = [celltmp.vLen;nEnd-nBeg];
                celltmp.vCen = [celltmp.vCen;nCount];
                nBeg = nEnd;
            end
            nEnd = nEnd + 1;
        else
            nCount = tc{i,j};
            celltmp.vCInt = [celltmp.vCInt;nBeg,nEnd];
            celltmp.nCInt = celltmp.nCInt + 1;
            celltmp.vLen = [celltmp.vLen;nEnd-nBeg];
            celltmp.vCen = [celltmp.vCen;nCount];
            nBeg = nEnd;
        end
    end
    if u <= 47
        if any(vRandB == u)
            dataPlacebo.nTrain = dataPlacebo.nTrain + 1;
            dataPlacebo.censor{dataPlacebo.nTrain} = celltmp;
        else
            dataPlacebo.nTest = dataPlacebo.nTest + 1;
            dataPlacebo.centest{dataPlacebo.nTest} = celltmp;
        end
    else
        if any(vRandA == u)
            dataDose.nTrain = dataDose.nTrain + 1;
            dataDose.censor{dataDose.nTrain} = celltmp;
        else
            dataDose.nTest = dataDose.nTest + 1;
            dataDose.centest{dataDose.nTest} = celltmp;
        end
    end
end
% save('dataBladder.mat','dataDose','dataPlacebo');
dataFull = cell(1,2);
dataFull{1} = ShrinkSet(dataDose);
dataFull{2} = ShrinkSet(dataPlacebo);
end

function dataFull = generateSkin()
% small function to generate Skin data set
t = readtable('skinCan.txt','ReadVariableNames',false);
[nrow,~] = size(t);
dataDose.nU = 143; % N1
dataDose.nTrain = 0;dataDose.nTest = 0;
dataDose.censor = cell(1,143);   % for training
dataDose.centest = cell(1,143);     % for testing
dataDose.XBeg = 0;dataDose.XEnd = 0;
dataDose2.nU = 143; % N2
dataDose2.nTrain = 0;dataDose2.nTest = 0;
dataDose2.censor = cell(1,143);   % for training
dataDose2.centest = cell(1,143);     % for testing
dataDose2.XBeg = 0;dataDose2.XEnd = 0;
vRandA = randperm(dataDose.nU);
vRandA = vRandA(1:round(dataDose.nU/2));

dataPlacebo.nU = nrow/3-143;  % N1
dataPlacebo.nTrain = 0;
dataPlacebo.nTest = 0;
dataPlacebo.censor = cell(1,dataPlacebo.nU);
dataPlacebo.centest = cell(1,dataPlacebo.nU);
dataPlacebo.XBeg = 0;dataPlacebo.XEnd = 0;
dataPlacebo2.nU = nrow/3-143; % N2
dataPlacebo2.nTrain = 0;
dataPlacebo2.nTest = 0;
dataPlacebo2.censor = cell(1,dataPlacebo2.nU);
dataPlacebo2.centest = cell(1,dataPlacebo2.nU);
dataPlacebo2.XBeg = 0;dataPlacebo2.XEnd = 0;
vRandB = dataDose.nU+randperm(dataPlacebo.nU);
vRandB = vRandB(1:round(dataPlacebo.nU/2));

for i = 1:nrow/3
    u = i;           % user id
    celltmp.vCInt = [];
    celltmp.nCInt = 0;
    celltmp.vLen = [];
    celltmp.vCen = [];
    celltmp.X = [];
    ctime = strsplit(t{3*i-2,1}{1},' ');ctime = ctime(2:end);vtime = arrayfun(@(x)str2double(x),ctime)/30; % time in days, too large
    cN1 = strsplit(t{3*i-1,1}{1},' ');cN1 = cN1(2:end);vN1 = arrayfun(@(x)str2double(x),cN1);
    cN2 = strsplit(t{3*i,1}{1},' ');cN2 = cN2(2:end);vN2 = arrayfun(@(x)str2double(x),cN2);
    
    celltmp.vCInt = [0,vtime(1:end-1);vtime]';
    celltmp.nCInt = length(vtime);
    celltmp.vLen = (vtime-[0,vtime(1:end-1)])';
    celltmp.vCen = vN1';
    celltmp.X = [];
    
    celltmp2.vCInt = [0,vtime(1:end-1);vtime]';
    celltmp2.nCInt = length(vtime);
    celltmp2.vLen = (vtime-[0,vtime(1:end-1)])';
    celltmp2.vCen = vN2';
    celltmp2.X = [];
    if length(celltmp.vCen) ~= length(celltmp.vLen) || length(celltmp2.vCen) ~= length(celltmp2.vLen)
        celltmp.vCen = celltmp.vCen(1:celltmp.nCInt);
        celltmp2.vCen = celltmp2.vCen(1:celltmp2.nCInt);
    end
    if u > 143
        dataPlacebo.XEnd = max(dataPlacebo.XEnd,max(vtime));
        dataPlacebo2.XEnd = dataPlacebo.XEnd;
        if any(vRandB == u)
            dataPlacebo.nTrain = dataPlacebo.nTrain + 1;
            dataPlacebo.censor{dataPlacebo.nTrain} = celltmp;
            dataPlacebo2.nTrain = dataPlacebo2.nTrain + 1;
            dataPlacebo2.censor{dataPlacebo2.nTrain} = celltmp2;
            
        else
            dataPlacebo.nTest = dataPlacebo.nTest + 1;
            dataPlacebo.centest{dataPlacebo.nTest} = celltmp;
            dataPlacebo2.nTest = dataPlacebo2.nTest + 1;
            dataPlacebo2.centest{dataPlacebo2.nTest} = celltmp2;
        end
    else
        dataDose.XEnd = max(dataDose.XEnd,max(vtime));
        dataDose2.XEnd = dataDose.XEnd;
        if any(vRandA == u)
            dataDose.nTrain = dataDose.nTrain + 1;
            dataDose.censor{dataDose.nTrain} = celltmp;
            dataDose2.nTrain = dataDose2.nTrain + 1;
            dataDose2.censor{dataDose2.nTrain} = celltmp2;
        else
            dataDose.nTest = dataDose.nTest + 1;
            dataDose.centest{dataDose.nTest} = celltmp;
            dataDose2.nTest = dataDose2.nTest + 1;
            dataDose2.centest{dataDose2.nTest} = celltmp2;
        end
    end
end
dataFull = cell(1,4);
dataFull{1} = ShrinkSet(dataDose);
dataFull{2} = ShrinkSet(dataDose2);
dataFull{3} = ShrinkSet(dataPlacebo);
dataFull{4} = ShrinkSet(dataPlacebo2);
end

function data = ShrinkSet(data)
% shrink the data set according to the actual size
data.censor = data.censor(1:data.nTrain);
data.centest = data.centest(1:data.nTest);
end

