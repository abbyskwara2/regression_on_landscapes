function plotSweepPanel(ax, glvData)
% muList sigmaList qual2
muList = glvData.muList;
sigmaList = glvData.sigmaList;
glv_qual2 = glvData.qual2;

%%
ax.Box = 'on';

im = imagesc(muList, sigmaList, glv_qual2, [0.5 1]);
colormap(parula(255));
im.AlphaData = ~isnan(glv_qual2);
ax.Color = 0.8*[1 1 1];
hold on;
xlabel('mu')
ylabel('sigma')
axis square
axis xy
adjustSizes(ax,1,13);
xlabel('\it\mu','FontSize',13);
ylabel('\it\sigma','FontSize',13);
title(sprintf('$\\gamma=%.1f, \\zeta=%.1f$', glvData.gamma, glvData.std_K),'FontSize',12,'FontWeight','normal','Interpreter','latex') % char(8199*[1 1]),

end