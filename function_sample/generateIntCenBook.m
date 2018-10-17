function data = generateIntCenBook(options)
%% Setting Censoring all of them
data.XBeg = 0;
data.XEnd = 60;
nLineWidth = 1.3;
if isfield(options,'f')
    f.x = linspace(data.XBeg,data.XEnd,options.testSample)';
    f.y = options.f;
    data.fBasic = @(x)linear(x,f);
    data.fileName = 'CenGP';
else
    switch options.dataType
        case 1
            data.fBasic = @(x)7*(x>=0 & x<10)+2*(x>=10 & x<20)+7*(x>=20 & x<30)...
                +2*(x>=30 & x<40)+7*(x>=40 & x<50)+2*(x>=50 & x<60);
            data.fileName = 'CenPL';
        case 2
            data.fBasic = @(x)10*exp(-x/15)+5*exp(-(x-35).^2/100);
            data.fileName = 'CenExpNormal';
        case 3
            data.fBasic = @(x)4*(1+sin(2*pi/20*x));
            data.fileName = 'CenSine';
        case 4
            data.fBasic = @(x)5*exp(-(x-30).^2/200);
            data.fileName = 'CenLinear';
        case 6
            data.fBasic = @(x) 100 * normpdf(x, 0.2, 0.05) + 80 * normpdf(x, 0.4, 0.03) + ...
                 60 * normpdf(x, 0.5, 0.03) + 40 * normpdf(x, 0.6, 0.03) + ...
                 20 * normpdf(x, 0.7, 0.03);
            data.XEnd = 1;
            data.fileName = 'CenNew';
        case 5 % 2-dimensional data
            data = generate2D(options);
            return
    end
end

vT = linspace(data.XBeg,data.XEnd,200)';
data.nU = options.censorExperiment/2;  % training number  
data.nTrain = data.nU;  % training number  
data.nTest = options.censorExperiment/2; % testing number
%% generate one sequence
nfMax = max(data.fBasic(vT))+2;
data.censor = cell(1,data.nU);      % censored train data set
data.train = cell(1,data.nU);       % complete train data set
data.test = cell(1,data.nU);       % complete test data set
data.centest = cell(1,data.nTest);  % censored test data set
%% equally censored all points 
for u = 1:options.censorExperiment
    ir = rand();
    if options.dataWeight ~= 0  % add variations
        fCurrent = @(x)(0.5+ir)*data.fBasic(x);
        nfMax = max(fCurrent(vT))+2;
        data.raw = sampleIPP(fCurrent,data.XBeg,data.XEnd,nfMax,1);
        nLenTotal = data.XEnd-data.XBeg;
    else                        % do not add variations
        fCurrent = data.fBasic;
        data.raw = sampleIPP(fCurrent,data.XBeg,data.XEnd,nfMax,1);
        nLenTotal = data.XEnd-data.XBeg;
    end

    data.dice = rand(size(data.raw))>0;
    % generate a censored ratio
    vRatio = gamrnd(ones(options.censorNum,1),1);     % ratio for generating the periods
    vRatio = vRatio./sum(vRatio);
    centmp.vCInt = zeros(options.censorNum,2);        % store the intervals
    centmp.nCInt = options.censorNum;                 % store the number of intervals
    centmp.vLen = zeros(options.censorNum,1);         % store the length of each interval
    centmp.vCen = zeros(options.censorNum,1);         % store the number of points in each interval
    centmp.fCurrent = fCurrent;
    centmp.nLike = 0;                                 % constant term in the likelihood
    centmp.X = []; % there are no un-censored points
    nBeg = 0;
    for ic = 1:options.censorNum
        nEnd = nBeg + nLenTotal*vRatio(ic);
        vFlag = data.raw >= nBeg & data.raw <= nEnd;  % can be faster here
        centmp.vCInt(ic,:) = [nBeg,nEnd];
        centmp.vLen(ic) = nLenTotal*vRatio(ic);
        nPoint = sum(vFlag);
        centmp.vCen(ic) = nPoint;
        nBeg = nBeg + nLenTotal*vRatio(ic);
    end
    centmp.nLike = computeConst(centmp.vCen);
    if u <= data.nU  % training set
        data.censor{u} = centmp;
        data.train{u}.X = data.raw';
    else             % testing set
        data.centest{u-data.nU} = centmp;
        data.test{u-data.nU}.X = data.raw';
    end
end
%% plotting
if options.dataPlot == 1
    fig = figure;
    set(fig,'Position',[0,0,700,600]);
    nPlot = min(data.nU,4);
    for u = 1:nPlot
        subplot(nPlot,1,u);hold on;
        for i = 1:data.censor{u}.nCInt
            rectangle('Position',[data.censor{u}.vCInt(i,1),-1,data.censor{u}.vLen(i),2],...
                'FaceColor',[0.75,0.75,1],'EdgeColor',[1,0,0],'LineWidth',nLineWidth);
        end
        for t = 1:length(data.train{u}.X)
            line(data.train{u}.X(t)*ones(1,2),[-1,1],'Color','k');
        end
        xlim([data.XBeg,data.XEnd]);% title(['Censored Sequence ',num2str(u)]);
        set(gca,'FontSize',13);
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        if u == data.nU
            xlabel('Time t');
        end
        
        ylabel(['Seq ',num2str(u-1)]);
    end
    % print figure
    set(gca,'XTick',unique([0,floor(data.censor{u}.vCInt(:,1))',60]));
    % title('Censor all');
    printFig(fig,'.\result_plot\CensorAll');
end
end

function vF = linear(vX,f)
% perform linear interpolation using f.x, f.y
vTx = f.x;
vTy = f.y;
if ~iscolumn(vTy)
    vTy = vTy';
end
nMin = min(vTx);
nMax = max(vTx);
nStep = vTx(2)-vTx(1);
vPos = (min(max(nMin,vX),nMax)-nMin)/nStep+1;
vPos_l = floor(vPos);
vPos_u = ceil(vPos);
vRatio_l = (vPos_u-vPos)./(vPos_u-vPos_l);
vRatio_u = (vPos-vPos_l)./(vPos_u-vPos_l);

vEqual = vPos_l == vPos_u;
vFequal = vTy(vPos_l);
vF = (vRatio_l.*vTy(vPos_l)+vRatio_u.*vTy(vPos_u));
vF(vEqual) = vFequal(vEqual);
vF(vX<nMin) = vTy(1);
vF(vX>nMax) = vTy(end);
end

function data = generate2D(options)
% for 2D data sets
switch options.dataType
    case 5
        data.XBeg = [0,0];
        data.XEnd = [60,60];
        data.fBasic = @(x,y)2*exp(-(x-40).^2/350-(x-40).*(y-40)/350-(y-40).^2/350);
end
data.nU = options.censorExperiment/5*4;
data.nTrain = data.nU;
data.nTest = options.censorExperiment/5;
data.train = cell(1,data.nTrain);
data.censor = cell(1,data.nTrain);
data.centest = cell(1,data.nTest);
for u = 1:data.nU
    fCurrent = data.fBasic;
    data.train{u}.X = sampleIPP2D(fCurrent,[data.XBeg(1),data.XEnd(1)],[data.XBeg(2),data.XEnd(2)]);
    data.censor{u}.fCurrent = fCurrent;
    % how to do censoring here
    
end
end