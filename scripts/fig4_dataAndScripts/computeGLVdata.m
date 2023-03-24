function glvData = computeGLVdata
fname = 'GLVdata.mat';
if exist(fname,'file')
    fprintf('Loading GLV simulation results from disk. (Delete to recompute from scratch.)\n');
    glvData = load(fname);
else
    TEST = false;
    params.S = 10;
    params.Teq = 10000;
    
    % Match params from Fig. S2 of Barbier et al, PNAS
    std_K = 0.3;
    abd_std = 0.1;
    N = params.S;
    
    params.growth = ones([N,1]);
    params.abd0   = exprnd(abd_std, [N,1]);
    params.propertyOfInterest = 'biomass';
    STOP_ON_FAILURE = true;
    
    if TEST
        repN = 2;
        muList = -1:2:3;
        sigmaList = 0.2:0.4:1;
        % muList = linspace(-1,3,9);
        % sigmaList = linspace(0.2,1,9);
    else    
        repN = 10;
        muList = linspace(-2,4,13);
        sigmaList = linspace(0.2,1.5,14);
    end
    [mumu, ss, rr] = meshgrid(muList, sigmaList, 1:repN);
    approxQual = NaN([N, size(mumu)]);
    biomass = NaN(size(mumu));
    landscape = NaN([2^N,size(mumu)]);
    
    tic
    for ii=1:length(mumu(:))
        fprintf('%d/%d\n',ii, length(mumu(:)))
        rng(ii);
        mu = mumu(ii);
        sigma =ss(ii);
        A = normrnd(mu/N, sigma/sqrt(N), [N,N]);
        A = tril(A,-1);
        params.inters = A+A';
        while true % get random K but ensure they are all positive
            params.K = 1 + std_K*randn([N,1]);
            if all(params.K>0), break, end
        end        
        GLV = cflLotkaVolterra(params);
    
        try
            divergenceDetected = false;
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
        end
    end
    toc
    meanApproxQual = mean(approxQual,4,'omitnan');
    qual1 = squeeze(meanApproxQual(1,:,:));
    qual2 = squeeze(meanApproxQual(2,:,:));
    
    divCount = squeeze(sum(any(isnan(approxQual)),4));
    repCount = size(approxQual,4);
    tooManyDivergences = divCount > repCount/2;
    qual2(tooManyDivergences)=NaN;
    save(fname, 'muList','sigmaList','qual2');
    glvData.muList = muList;
    glvData.sigmaList = sigmaList;
    glvData.qual2 = qual2;
    
end % if exist fname
