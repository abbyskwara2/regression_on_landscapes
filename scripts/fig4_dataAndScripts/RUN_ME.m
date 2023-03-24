function RUN_ME
% Generate the simulation figure for Skwara et al.
% Simulation code: Mikhail Tikhonov

    % Simulation data for panel A
    glvData = computeGLVdata;
    bdryCurve = phaseBoundaryCurve(Inf);
    % Simulation data for panel B
    crmData = computeCRMdata;

    plotFigure(glvData, bdryCurve, crmData);
    saveas(gcf,'simFig.fig');
    saveas(gcf,'simFig.pdf');
end

