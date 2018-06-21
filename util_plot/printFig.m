function printFig(fig,strFile)
% Print the figure into PDF file
% Input :    fig     -- figure
%           strFile -- output file position
% Output:   None


set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(fig,strFile,'-dpdf','-r0');
end