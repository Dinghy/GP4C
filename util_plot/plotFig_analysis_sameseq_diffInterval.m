clear;close all;clc;
addpath(genpath('./'));
%%
% obtain a sample from a point process
XBeg = 0;
XEnd = 60;
nlen = XEnd-XBeg;
fBasic = @(x)5*ones(size(x));
vt = linspace(XBeg,XEnd,200)';

vnum = 1:10;
nrepeat = 300;
mlike = zeros(length(vnum),nrepeat);
vdata = sampleIPP(fBasic,XBeg,XEnd,max(fBasic(vt))+2,1);
vlogcache = zeros(1,length(vdata)+1);

% ln(m!), using vlogcache(i) to save ln(i!),i=0,...,m
for i = 2:length(vdata)+1
    vlogcache(i) = vlogcache(i-1)+log(i-1);
end

for n = vnum
    for j =1:nrepeat        
        vRatio = gamrnd(ones(n,1),1);
        vRatio = vRatio./sum(vRatio);
        nBeg = 0;
        nlike = 0;
        vnum = zeros(1,n);  % record the number of points in the censored data
        for k = 1:n
            nEnd = nBeg + nlen*vRatio(k);
            nnum = sum(vdata <= nEnd & vdata > nBeg);
            nlike = nlike + nnum*log(5*(nEnd-nBeg));
            vnum(k) = nnum;
            nBeg = nEnd;
        end
        mlike(n,j) = nlike-5*nlen-computeConst(vnum);
    end
end
%%
fig = figure;hold on;grid on;box on;
boxplot(mlike');
% vx = 1:10;
% vy = length(vdata)*log(nlen*5./vx)-vx.*gammaln(1+length(vdata)./vx)-300;
% plot(vx,vy,'r*-');
set(gca,'FontSize',13);
xlabel('The number of the intervals N');
ylabel('Logarithm of the likelihood');
printFig(fig,'result_plot\LogLikelihood'); % Figure 11