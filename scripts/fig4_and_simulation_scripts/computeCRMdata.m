function crmData = computeCRMdata
fname = 'CRMdata.mat';
if exist(fname,'file')
    fprintf('Loading CRM simulation results from disk. (Delete to recompute from scratch.)\n')
    crmData = load(fname);
    crmData.L = 8;
else
    TEST = false;
    LEAK_TO_NEXT_LAYER = 0.5;
    dLEAK = 0.05;
    L = 8;
    S = 10;
    if TEST 
        repN = 2;
        listLambda = linspace(1/L, 1, 3);
    else
        repN = 20;
        listLambda = linspace(1/L, 1, 15);
    end
    params.equilibrationTime = 1e4;
    params.m = ones(S,1);
    params.propertyOfInterest = 'lastResource';
    %%
    [lamlam, rr] = meshgrid(listLambda, 1:repN);
    approxQual = NaN([S, size(rr)]);
    landscape = NaN([2^S,size(rr)]);
    totalIter = length(rr(:));
    minlen = NaN(size(lamlam));
    
    for ii=1:totalIter
        fprintf('%d/%d\n',ii,totalIter)
        seed = ii;
        rng(seed);
        lambda = lamlam(ii);
    
        params.K = [1e6, zeros(1,L-1)]; % only first resource is supplied externally
    
        [circuit, iter] = crossFeedingCircuitInstance(S,L,lambda,seed);
        minlen(ii) = length(circuit.G.shortestpath(1,L))-1;
    
        params.D = circuit.D*LEAK_TO_NEXT_LAYER.*(1+dLEAK*randn(size(circuit.D)));
        params.C = circuit.C;
        % compute the landscape
        CRM = cflCRM_crossFeeding(params);
        f = CRM.computeAll;
        landscape(:,ii) = f;
        c = landscape2fourier(f);
        p = fourier2power(c);
        p = p(2:end);
        p = p/sum(p);
        approxQual(:,ii) = cumsum(p);
    end
    save(fname,'approxQual','minlen');
    crmData.approxQual = approxQual;
    crmData.minlen = minlen;
    crmData.L = L;
end % if exist fname
end
