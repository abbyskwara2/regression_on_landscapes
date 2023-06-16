function RUN_ME
% Generate the simulation figure for Skwara et al.
% Simulation code by Mikhail Tikhonov <tikhonov@wustl.edu>
%
% Figure plotting uses "arrow.m" by Erik A. Johnson <JohnsonE@usc.edu>
% http://www.usc.edu/civil_eng/johnsone/


    % Simulation data for Fig 4A
    bdryCurve = phaseBoundaryCurve(Inf);
    % Simulation data for Fig 4B
    crmData = computeCRMdata;

    % Match params from Fig. S2 of Barbier et al, PNAS
    stdK = 0.3; gamma = 1;

    params.PLOT_DIAGNOSTIC = true;
    params.COARSE = false;
    glvData = computeGLVdata(stdK,gamma, params);
    plotFigure(glvData, bdryCurve, crmData);

    saveas(gcf,sprintf('simFig_stdK=%g_gamma=%g.fig', stdK, gamma));
    saveas(gcf,sprintf('simFig_stdK=%g_gamma=%g.pdf', stdK, gamma));

    %% Dig into a pixel: distribution of scores across many replicates?
    repN = 1000;
    mu = 1;
    sigma = 0.8;
    manyReplicatesData = characterizeDistributionAcrossReplicates(mu,sigma,gamma,stdK, repN);
    plotAtypicalExamples(glvData,bdryCurve,manyReplicatesData);
    saveas(gcf,'simFig_rugged.fig');
    saveas(gcf,'simFig_rugged.pdf');

    %% make SI figure with 4d parameter sweep    
    % Compute & save diagnostics
    stdKlist = [0.1 0.3 1];
    gammaList = -1:0.5:1;

    params.PLOT_DIAGNOSTIC = true;
    params.COARSE = false;
    % Sweep parameters     
    for stdK = stdKlist
        for gamma = gammaList
            computeGLVdata(stdK,gamma, params);
        end
    end

    %% Plot the actual figure
    W = 33;
    H = 20;
    clf;
    set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters',...
        'PaperSize', [W H], 'PaperPosition',[0 0 W H],...
        'Units','Centimeters','Position',[2 2 W H]); 
    params.PLOT_DIAGNOSTIC = false;
    tiledlayout(length(stdKlist),length(gammaList),"TileSpacing","compact");
    for stdK = stdKlist
        for gamma = gammaList
            glvData = computeGLVdata(stdK,gamma, params);
            plotSweepPanel(nexttile, glvData);
        end
    end
    cbh = colorbar(gca); 
    % To position the colorbar as a global colorbar representing all tiles 
    cbh.Layout.Tile = 'east';
    saveas(gcf,'simFig_SI.fig');
    saveas(gcf,'simFig_SI.pdf');
    %%
end



