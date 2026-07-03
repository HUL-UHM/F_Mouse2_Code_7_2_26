% setfigparms - script to setup the correct figure parameters
% assumes you assigned 'colordef none' before opening figures

set(gcf, 'Color', [1 1 1]); 

h=gca;
set(get(h,'YLabel'),'FontSize',24)%Used to be 25
set(get(h,'XLabel'),'FontSize',24)
set(get(h,'Title'),'FontSize',28)
%set(get(h,'Title'),'Color',[1 1 1])
%set(h,'Position',[0.19 0.24 .6 .6]); % This provides a slide mount margin
%set(gca, 'Color', [0 0 .6])
%set(gca, 'Color', [1 1 1])
set(gca, 'FontSize', 20) % used to be 16
set(gca, 'LineWidth', 4)
hh=get(gca,'Children');
%set(hh, 'Color', [1 1 0])
%set(hh, 'LineWidth' ,2) 
%set(gca, 'XColor', [1 1 1])
%set(gca, 'YColor', [1 1 1])
%set(gcf, 'Color', [0 0 .6])
set(gcf,'invertHardcopy', 'off')
%disp('print  -dtiffnocompression  fname')

%for white background
%set(gcf,'Color',[1 1 1]) 
%to remove axis marks
%set(gca, 'XTickLabel', [])
%set(gca, 'YTickLabel', [])
