function plotAtypicalExamples(glvData,bdryCurve,manyReplicatesData)
% muList sigmaList qual2
muList = glvData.muList;
sigmaList = glvData.sigmaList;
glv_qual2 = glvData.qual2;

%%
W = 22;
H = 7.6;
w0 = 1.6;
h0 = 1.4;
wh = 5;
dw = 1.5;
overlayCol = 'k';
lw = 0.5;

clf;
set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters','Position',[4 4 W H]); 
axA = axes('Units','Centimeters','Position',[w0, h0, wh+1.5, wh]); 
box on

im = imagesc(muList, sigmaList, glv_qual2, [0.5 1]);
colormap(parula(255));
im.AlphaData = ~isnan(glv_qual2);
axA.Color = 0.8*[1 1 1];
hold on;
plot(bdryCurve(1,:), bdryCurve(2,:), [overlayCol,'--'],'LineWidth',lw);
xlabel('mu')
ylabel('sigma')
axis square
axis xy
colorbar
adjustSizes(axA,1,12);
%xlabel('Interaction strength $\mu$','Interpreter','LaTeX','FontSize',16);
xlabel('Strength of interactions \it(\mu)','FontSize',13);
ylabel('Std of interactions \it(\sigma)','FontSize',13);
%ylabel('$\sigma$','Interpreter','LaTeX','FontSize',16);
title('     Mean R^2 of 2nd-order model','FontSize',12,'FontWeight','normal') % char(8199*[1 1]),

% The test asterisk looks better than Matlab's marker. But text position is
% not precisely centered where requested. Hence the fudge parameter -0.035.
% Uncomment next line to check that the alignment is indeed correct.
%plot(manyReplicatesData.mu, manyReplicatesData.sigma, 'k*', 'MarkerSize',20);
text(manyReplicatesData.mu, manyReplicatesData.sigma-0.035, '*', 'FontSize',26,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

w0 = w0+wh+1.5+dw;
axB = axes('Units','Centimeters','Position',[w0, h0, wh, wh]); 
box on
hold all
axis square
q2 = manyReplicatesData.approxQual(:,2);
h = histogram(q2,0.602:0.005:1,'FaceColor','k');
[q2srt, order] = sort(q2);
axB.YScale = 'log';
axB.YLim = [0.8 300];
axB.YTick = [1 10 100 1000];
axB.XLim = [0.6 1.005];
hold on
%histogram(q2srt(1:3), h.BinEdges);

for ii=1:2
    arrow([q2srt(ii), 1.2], [q2srt(ii), 1.1],'Color','r');
end
text(axB.XLim(1)+0.04, 2.3, {'rugged','landscapes'},'HorizontalAlignment','left','FontSize',12, 'VerticalAlignment','bottom');
labelX = mean(axB.XLim)+0.14;
labelY = 80;
text(labelX, labelY, {'smooth','landscapes'},'HorizontalAlignment','right','FontSize',12, 'VerticalAlignment','middle');
%arrow([labelX-0.01, labelY], [labelX, labelY]);

title('Distribution across replicates','FontSize',12,'FontWeight','normal')
xlabel('R^2 of second-order model')
adjustSizes(axB,1,12);

w0 = w0+wh+dw+0.5;
axC = axes('Units','Centimeters','Position',[w0, h0, wh, wh]); 
box on
hold all
axis square
adjustSizes(axC,1,12);
qual = manyReplicatesData.approxQual;
powerDistrib = diff([zeros(size(qual,1),1), qual],[],2);
h1 = plotPatch(mean(powerDistrib,"omitmissing"), std(powerDistrib,"omitmissing"), 'k');
%plot(mean(powerDistrib,"omitmissing"),'k-','LineWidth',1);
hold on
h2 = plot(powerDistrib(order(1:2),:)','r-','LineWidth',1);
axC.XLim = [1 5];  
axC.YLim = [0 0.8];  
title('Fraction variance explained','FontWeight','normal');
lg = legend([h1, h2(1)], 'Typical landscapes','Highlighted in B','FontSize',10);
lg.Position(1) = lg.Position(1)+0.008;
lg.Box = 'off';


xlabel('Order')
labelYOffset = 0.1;
text(axA, -1.1,wh+labelYOffset,'A','Units','Centimeters','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','bottom')
text(axB, -1.1,wh+labelYOffset,'B','Units','Centimeters','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','bottom')
text(axC, -1.1,wh+labelYOffset,'C','Units','Centimeters','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','bottom')

set(gcf, 'Color','w');
end

function h = plotPatch(mn, er, col)
    xs = 1:length(mn);
    h = plot(xs,mn,'-','Color', col, 'LineWidth',2);
    patch([xs, fliplr(xs)],[mn+er, fliplr(mn-er)],col,'FaceAlpha',0.5,'EdgeColor','none');
end
