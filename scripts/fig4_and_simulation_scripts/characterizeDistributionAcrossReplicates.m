function data = characterizeDistributionAcrossReplicates(mu,sigma,gamma,std_K, repN)
fname = sprintf('manyReplicates_mu=%g_sigma=%g_gamma=%g_stdK=%g', mu, sigma, gamma, std_K);
if exist([fname, '.mat'],'file')
    fprintf(sprintf('Loading many-replicate analysis from disk. (Delete %s.mat to recompute from scratch.)\\n',fname));
    data = load([fname, '.mat']);
    data = data.data;
else

    N = 10;
    abd0_std = 0.1;
    STOP_ON_FAILURE = true;
    
    approxQual = NaN([repN, N]);

    tic
    for rr=1:repN
        fprintf('%d/%d\n',rr, repN)
        rng(rr);
        
        sigma_eta = sqrt((1+gamma)/2)*sigma/sqrt(N);
        sigma_xi = sqrt((1-gamma)/2)*sigma/sqrt(N);
        eta = tril(normrnd(0, sigma_eta, [N,N]),-1);
        xi = tril(normrnd(0, sigma_xi, [N,N]),-1);
        aij = eta+xi;
        aji = (eta-xi)';
        params = [];
        params.Teq = 10000;    
        params.S = N;    
        params.growth = ones([N,1]);
        params.abd0   = exprnd(abd0_std, [N,1]);
        params.propertyOfInterest = 'biomass';
        
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
            f = GLV.computeAll(STOP_ON_FAILURE);
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
            approxQual(rr,:) = cumsum(p);
        end
    end
    toc
    
    data.approxQual = approxQual;
    data.mu = mu;
    data.sigma = sigma;
    data.std_K = std_K;
    data.gamma = gamma;
    save([fname, '.mat'], 'data');
    
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