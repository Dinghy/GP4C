function plotEbar(vX,mY,fig,vColor,vPos,options)
% Plot an error bar on the figure
% Input :  vX   -- plot x, a vector
%          mY   -- plot y, a matrix with the column number equals to
%                  the element number in vX
%          fig  -- figure
%          options  -- errorbar, marker
%          v.c. -- color of the plot, a vector with length 3
%          v.p. -- position of the plot, a vector with length 4
% Output:  None
nWidth = 1.5;
figure(fig);
subplot('Position',vPos);
if ~isvector(vX)
    error('Illegal inputs for vX');
end
if ~isvector(vColor) || length(vColor) ~= 3
    error('Illegal color');
end

if isvector(mY)
    nMeanY = mean(mY);
    nSigY = sqrt(cov(mY));
    hold on;
    if options.errbar == 1
        nMed = median(mY);
        nLow = nMed-quantile(mY,0.25);
        nHigh = quantile(mY,0.75)-nMed;
        [nMed,nLow,nHigh]
        rectangle('Position',[min(vX),nMed-nLow,max(vX)-min(vX),nHigh+nLow],'FaceColor',(vColor+[1,1,1])/2,'EdgeColor','w');
        plot(vX,nMeanY*ones(size(vX)),'Color',vColor,'LineWidth',nWidth);
    else
        if options.marker == 1
            plot(vX,nMeanY*ones(size(vX)),'Marker','*','Color',vColor,'LineWidth',nWidth);
        else
            plot(vX,nMeanY*ones(size(vX)),'Color',vColor,'LineWidth',nWidth);
        end
    end
else
    vMeanY = mean(mY,1);
    hold on;
    if options.errbar == 1
        vMedY = median(mY,1);
        vHighY = quantile(mY,0.75,1)-vMedY;
        vLowY = vMedY - quantile(mY,0.25,1);
        errorbar(vX,vMedY,vLowY,vHighY,'Color',vColor,'LineWidth',nWidth);
    else
        if options.marker == 1
            plot(vX,vMeanY,'Marker','*','Color',vColor,'LineWidth',nWidth);
        else
            plot(vX,vMeanY,'Color',vColor,'LineWidth',nWidth);
        end
    end
end