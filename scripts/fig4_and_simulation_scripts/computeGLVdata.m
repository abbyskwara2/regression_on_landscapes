function glvData = computeGLVdata(std_K, gamma, params)
if nargin<3
    params = [];
end
[PLOT_DIAGNOSTIC, COARSE] = parseParams(params);
fname = sprintf('GLVdata_stdK=%g_gamma=%g', std_K, gamma);
if COARSE
    fname = [fname, '_COARSE'];
end
if exist([fname, '.mat'],'file')
    fprintf('Loading GLV simulation results from disk. (Delete to recompute from scratch.)\n');
    glvData = load([fname, '.mat']);
    glvData = glvData.glvData;
else

    params.S = 10;
    params.Teq = 10000;
    
    abd0_std = 0.1;
    N = params.S;
    
    params.growth = ones([N,1]);
    params.abd0   = exprnd(abd0_std, [N,1]);
    params.propertyOfInterest = 'biomass';
    STOP_ON_FAILURE = true;
    
    if COARSE
        repN = 3;
        muList = linspace(-1,3,9);
        sigmaList = linspace(0.2,1,9);
    else    
        repN = 10;
        muList = linspace(-2,4,13);
        sigmaList = linspace(0.2,1.5,14);
    end
    [mumu, ss, rr] = meshgrid(muList, sigmaList, 1:repN);
    approxQual = NaN([N, size(mumu)]);
    biomass = NaN(size(mumu));
    landscape = NaN([2^N,size(mumu)]);

    % to double-check the code: compute the 4 statistics of the LV system
    % that Barbier et al uses as control parameters. These are the
    % parameters that we sweep (so the LV matrix is *supposed* to have the
    % correct stats); but it never hurts to independently check!
    parameterVerification = NaN([4, size(mumu)]);
    spotCheck = NaN(size(mumu)); % independent check of Fourier
    explainedVariance2 = NaN(size(mumu)); % independent check of Fourier
    tic
    for ii=1:length(mumu(:))
        fprintf('%d/%d\n',ii, length(mumu(:)))
        mu = mumu(ii);
        sigma =ss(ii);
        rep = rr(ii);
        seed = round(rep + sigma*1e4 + (mu+1)*1e6 + std_K*1e8 + (gamma+1)*1e9);
        rng(seed);
        
        sigma_eta = sqrt((1+gamma)/2)*sigma/sqrt(N);
        sigma_xi = sqrt((1-gamma)/2)*sigma/sqrt(N);
        eta = tril(normrnd(0, sigma_eta, [N,N]),-1);
        xi = tril(normrnd(0, sigma_xi, [N,N]),-1);
        aij = eta+xi;
        aji = (eta-xi)';
        params.inters = mu/N + aij + aji;
        
        while true % get random K but ensure they are all positive
            %params.K = 1 + std_K*randn([N,1]);
            % gamma(a,b) has mean ab and variance a*b^2
            % so to have mean 1, b=1/a.
            % and to have variance stdK^2, a*(1/a)^2=stdK^2 -> a=1/stdK^2
            a = 1/std_K^2;
            params.K = gamrnd(a,1/a,[N,1]);
            if all(params.K>0), break, end
        end        
        GLV = cflLotkaVolterra(params);
    
        try
            divergenceDetected = false;
            parameterVerification(:,ii) = GLV.stats();
            f = GLV.computeAll(STOP_ON_FAILURE);
            landscape(:,ii) = f;
            biomass(ii) = GLV.getF(true(1,params.S));
        catch ME
            if strcmpi(ME.identifier, 'getF:fail')
                divergenceDetected = true;
                fprintf('Divergence detected\n');
            else
                rethrow(ME);
            end
        end
    
        if ~divergenceDetected
            c = landscape2fourier(f);
            p = fourier2power(c);
            % ignore the constant term; only look at % variance explained
            p = p(2:end); 
            p = p/sum(p);
            approxQual(:,ii) = cumsum(p);

            [spotCheck(ii), explainedVariance2(ii)] = runSpotChecks(f,c);
        end
    end
    toc
    fprintf('Independently verificaiton of Fourier coefficient calculation.\n\tMismatch (shuold be close to 0): %g\n', ...
        max(abs(spotCheck(:))));
    explainedVariance2_also = squeeze(approxQual(2,:,:,:));
    fprintf('Independently verificaiton of predictive power of 2nd order model.\n\tMismatch (should be close to 0): %g\n', ...
        max(abs(explainedVariance2(:)-explainedVariance2_also(:))));

    meanApproxQual = mean(approxQual,4,'omitnan');
    qual1 = squeeze(meanApproxQual(1,:,:));
    qual2 = squeeze(meanApproxQual(2,:,:));
    
    divCount = squeeze(sum(any(isnan(approxQual)),4));
    repCount = size(approxQual,4);
    if COARSE
        tooManyDivergences = divCount > 0;
    else
        tooManyDivergences = divCount > repCount/2;
    end
    qual2(tooManyDivergences)=NaN;
    glvData.approxQual = approxQual;
    glvData.muList = muList;
    glvData.sigmaList = sigmaList;
    glvData.qual2 = qual2;
    glvData.std_K = std_K;
    glvData.gamma = gamma;
    save([fname, '.mat'], 'glvData');
    
    if PLOT_DIAGNOSTIC
        %% save diagnostic plots
        f = figure;
        makeDiagnosticPlots(parameterVerification);
        print([fname, '_paramsVerification.png'], '-dpng','-r300');
        close(f);
    end
end % if exist fname
end

function makeDiagnosticPlots(parameterVerification)
    W=25;
    H=12;
    set(gcf, 'PaperPositionMode','Manual', 'PaperUnits','Centimeters','PaperSize', [W H], 'PaperPosition',[0 0 W H],'Units','Centimeters', 'Position',[0 1 W H]); 

    tiledlayout(2,4);
    paramNames = {'\mu','\sigma','\gamma', 'std K'};
    ranges = {[-2 4], [0 1.5], [-1 1], [0 1]};
    for pp = 1:4
        nexttile(pp); 
        colormap parula
        imagesc(squeeze(mean(parameterVerification(pp,:,:,:),4)), ranges{pp});
        title(paramNames{pp}); 
        colorbar
        axis xy
        axis square

        nexttile(pp+4);
        colormap parula
        imagesc(squeeze(std(parameterVerification(pp,:,:,:),[],4)), [0 1]);
        title(['std of ',paramNames{pp}]); 
        colorbar
        axis xy
        axis square
    end
end

function [spotCheck, explainedVariance2] = runSpotChecks(f,c)
    % spot-check a random Fourier coefficient
    S = log2(length(f));
    cube = fullCube(S);
    pm1 = 2*single(cube)-1;
    check = randi(length(c)); 
    split = prod(pm1(:, cube(check,:)),2); % partitions into + and -
    expected = (sum(f(split==1))-sum(f(split==-1)))/2^S;
    spotCheck = c(check) - expected;
    % explicitly check the explanatory power of second order model
    modelC = c;
    modelC(sum(cube,2)>2) =  0;
    modelF = fourier2landscape(modelC);
    explainedVariance2 = 1-var(f-modelF)/var(f);
end

function [PLOT_DIAGNOSTIC, COARSE] = parseParams(params)
    if isfield(params,'PLOT_DIAGNOSTIC')
        PLOT_DIAGNOSTIC = params.PLOT_DIAGNOSTIC;
    else
        PLOT_DIAGNOSTIC = true;
    end
    if isfield(params,'COARSE')
        COARSE = params.COARSE;
    else
        COARSE = false;
    end
end