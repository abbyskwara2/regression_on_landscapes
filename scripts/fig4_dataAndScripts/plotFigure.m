function plotFigure(glvData, bdryCurve, crmData)
% muList sigmaList qual2
muList = glvData.muList;
sigmaList = glvData.sigmaList;
glv_qual2 = glvData.qual2;

L = crmData.L;
crm_approxQual = crmData.approxQual;
crm_minlen = crmData.minlen;
crm_qual1 = squeeze(crm_approxQual(1,:,:));
crm_qual2 = squeeze(crm_approxQual(2,:,:));
crm_qual3 = squeeze(crm_approxQual(3,:,:));
crm_mnQual1 = accumarray(crm_minlen(:), crm_qual1(:),[L-1,1],@mean)';
crm_mnQual2 = accumarray(crm_minlen(:), crm_qual2(:),[L-1,1],@mean)';
crm_mnQual3 = accumarray(crm_minlen(:), crm_qual3(:),[L-1,1],@mean)';
crm_stdQual1 = accumarray(crm_minlen(:), crm_qual1(:),[L-1,1],@std)';
crm_stdQual2 = accumarray(crm_minlen(:), crm_qual2(:),[L-1,1],@std)';
crm_stdQual3 = accumarray(crm_minlen(:), crm_qual3(:),[L-1,1],@std)';

%%
W = 19;
H = 7.6;
w0 = 1.6;
h0 = 1.4;
wh = 5;
dw = 2;
overlayCol = 'k';
lw = 2;

clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[4 4 W H]); 
axA = axes('Units','Centimeters','Position',[w0, h0, wh+1.5, wh]); 
box on

im = imagesc(muList, sigmaList, glv_qual2, [0.5 1]);
colormap(parula(255));
im.AlphaData = ~isnan(glv_qual2);
axA.Color = 0.8*[1 1 1];
hold on;
plot(bdryCurve(1,:), bdryCurve(2,:), [overlayCol,'-'],'LineWidth',lw);
xlabel('mu')
ylabel('sigma')
axis square
axis xy
colorbar
adjustSizes(axA,1,13);
%xlabel('Interaction strength $\mu$','Interpreter','LaTeX','FontSize',16);
xlabel('Strength of interactions \it(\mu)','FontSize',13);
ylabel('Std of interactions \it(\sigma)','FontSize',13);
%ylabel('$\sigma$','Interpreter','LaTeX','FontSize',16);
title('Predicting total biomass','FontSize',12,'FontWeight','normal') % char(8199*[1 1]),

col = lines(3);

extraWA = 3;
w0 = w0+wh+dw;
w0 = w0+dw+extraWA;
axC = axes('Units','Centimeters','Position',[w0, h0, wh, wh]); 
box on
hold all
h1 = plotPatch(crm_mnQual1, crm_stdQual1, col(1,:));
h2 = plotPatch(crm_mnQual2, crm_stdQual2, col(2,:));
h3 = plotPatch(crm_mnQual3, crm_stdQual3, col(3,:));
axis([1 L-1 0 1]);
axis square
xlabel('Shortest chain length (\itL)');
ylabel('Variance explained (\itR ^2\rm)');
adjustSizes(axC,1,13);
axC.YTick = [0 1];
title('Predicting metabolite flux','FontSize',12,'FontWeight','normal')

[lg,icons] = legend([h1 h2 h3],{'','',''},'FontSize',10,'Orientation','horizontal');
lg.Position(2) = lg.Position(2)-0.02;
lg.Box = 'off';
wid = 0.07;
dwid = 0.05;
icons(4).XData = 1+[-wid 0]-2*(wid+dwid);
icons(6).XData = 1+[-wid 0]-(wid+dwid);
icons(8).XData = 1+[-wid 0];
text(axC, wh-0.05, wh-0.04,'Model order ','Units','Centimeters','FontSize',10,'VerticalAlignment','top',...
    'HorizontalAlignment','right');
text(axC, wh-0.1, wh-0.7,'1^{st}  2^{nd}  3^{rd}','Units','Centimeters','FontSize',10,'VerticalAlignment','top',...
    'HorizontalAlignment','right');


text(axA, -1.1,wh+0.5,'A','Units','Centimeters','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','bottom')
text(axA, 0.2,wh-0.2,'2^{nd} order','Units','Centimeters','FontSize',10,'VerticalAlignment','top')
text(axA, wh+1.8,wh+0.5,'B','Units','Centimeters','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','bottom')
text(axA, 7.35, wh+0.9,'Resource competition model with cross-feeding','Units','Centimeters','FontSize',12,'FontWeight','bold');
text(axA, 0, wh+0.9,'Generalized Lotka-Volterra','Units','Centimeters','FontSize',12,'FontWeight','bold');
text(axA, wh+0.2, wh+0.05,'\it R^2','Units','Centimeters','FontSize',12,'FontWeight','normal','VerticalAlignment','bottom');

set(gcf,'InvertHardCopy','Off');
set(gcf, 'Color','w');
end

function h = plotPatch(mn, er, col)
    xs = 1:length(mn);
    h = plot(xs,mn,'-','Color', col, 'LineWidth',2);
    patch([xs, fliplr(xs)],[mn+er, fliplr(mn-er)],col,'FaceAlpha',0.5,'EdgeColor','none');
end
