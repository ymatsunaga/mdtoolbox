function formatplot()
%% formatplot
% fomart handle graphics (fonts, lines, etc.) of the current figure
%
%% Syntax
%# formatplot;
%
%% Description
% Changes the handle graphics (fonts, lines, etc.) of the current
% figure to appropriate sizes, widths, colors.
%
%% Example
%# plot(randn(100,1), '-');
%# formatplot;
%# xlabel('step'); ylabel('value');
%
%% See also
% exportas
%

%% Root Object

%% Figure Object
set(gcf,'Color','w')
% set(gcf,'Resize','on')
% set(gcf,'PaperType','A4')
% set(gcf,'PaperOrientation','landscape')
% set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperUnits','normalized')
% %set(gcf,'PaperPosition',[0.2 0.2 0.6 0.6])
% set(gcf,'PaperPosition',[0 0 0.6 0.6])

%% Axes Object
set(gca,'GridLineStyle',':')
set(gca,'LineWidth',1)
set(gca,'FontSize',25)
% set(gca,'LineStyleOrder','k-o|k-x|k-+|k-*|k-s|k-d|k-v|k-^|k-<|k->|k-p|k-h')
set(gca,'LineStyleOrder','-|-o|-x|-+|-*|-s|-d|-v|-^|-<|->|-p|-h')
% set(gca,'LineWidth',2)
%set(gca,'FontName','Helvetica')
% h = get(gca,'children');
% set(h,'LineStyleOrder','k-o|k-x|k-+|k-*|k-s|k-d|k-v|k-^|k-<|k->|k-p|k-h')

%% Line Object
h=get(gca,'children');
% set(h,'MarkerEdgeColor','k')
%set(h,'MarkerFaceColor','w')
set(h,'MarkerSize',8)
% set(h,'FontSize',18)
% set(gca,'XTick',[10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2),10^(3)])
set(h,'LineWidth',2)

