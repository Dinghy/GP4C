%% Plot of Figure 3 in the paper
% Influence of b in Lemma 2.
clear;close all;clc;
load('gzdgz.mat');
addpath(genpath('./'));
%% left plot
vX = 10.^linspace(-3,3,5000)';
vX = vX.^2;
vYTrue = - queryGz(vX/2,gz,dgz);

% vYL1 = log(vX.^2+1);% vYL3 = log(vX.^2-0.5)+log(2)-0.5772156649;
vYL1 = log(vX+1);
vYL2 = log(vX+0.30);

vb = linspace(0,1,50);
vcov = zeros(size(vb));
vdX = 10.^linspace(-3,3,5000)';
vdX = vdX.^2;
vdYTrue = - queryGz(vdX/2,gz,dgz);
for i = 1:length(vb)
    vdYt = log(vdX+vb(i));
    vcov(i) = cov(vdYTrue-vdYt);
end
[~,iMin] = min(vcov);

fig = figure;set(fig,'Position',[100,100,600,350]);

subplot('Position',[0.08,0.16,0.4,0.8]);box on;
semilogx(vX,vYTrue,'r','LineWidth',1.5);hold on;
semilogx(vX,vYL1,'b','LineWidth',1.5);
semilogx(vX,vYL2,'m','LineWidth',1.5);
leg = legend('-G(-$$\varphi$$/2)','ln($$\varphi$$+1)','ln($$\varphi$$+0.3)','Location','best');
set(leg,'Interpreter','latex')
xlim([0.01,100]);grid on;
set(gca,'XTick',[0.01,0.1,1,10,100]);
ylabel('Function Value'); xl = xlabel('Value of $$\varphi$$');set(xl,'Interpreter','latex');set(gca,'FontSize',13);

subplot('Position',[0.59,0.16,0.4,0.8]);box on;
semilogy(vb,vcov,'LineWidth',1.5);hold on;
semilogy(vb(iMin),vcov(iMin),'ro','LineWidth',1.5);
ylabel('Empirical Variance');xlabel('Value of b');
set(gca,'XTick',unique([0,0.1,0.5:0.2:1,vb(iMin)]));
grid on;set(gca,'FontSize',13);
printFig(fig,'result_plot\HyperBound.pdf')


