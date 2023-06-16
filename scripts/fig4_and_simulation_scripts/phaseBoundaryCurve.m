% use Bunin PRE (2017) to compute the phase transition line
function xy = phaseBoundaryCurve(S,gamma)
% denote Z(mu, sigma) = sigma*h+mu (using same notation as Bunin 2017)
% we are looking for zeros of Z
if nargin<2
    gamma = 1;
end
if nargin<1
    S = Inf;
end
u = @(mu,sigma) (1-mu/S)./sigma;


deltas = linspace(-6,6,1000);
us = uOfDelta(deltas, gamma);
Delta = @(u)interp1(us,deltas,u);
% %% Sanity check: should be zero
% max(abs(W2(deltas)-deltas.*W1(deltas)-W0(deltas)))

%%
h = @(d)d.*sqrt(W2(d))./W1(d);

Z = @(mu,sigma) mu + sigma.*h(Delta(u(mu,sigma)));

%%
muVec = linspace(-2,4.5,100);
sigmaVec = linspace(0.1,1.4,100);
[mumu, sigsig] = meshgrid(muVec,sigmaVec);

zs = reshape(Z(mumu(:), sigsig(:)),size(mumu));
imagesc(muVec, sigmaVec, zs, [-1 1]);
axis xy

ms = [sigmaVec(:),sigmaVec(:)];
hold all
cc = contour(mumu, sigsig, zs,[0,0]);
xy = cc(:,2:end);
plot(xy(1,:), xy(2,:),'k-')
end


function u = uOfDelta(d, gamma)
% w0 = W0(d);
% sqrtw2 = sqrt(W2(d));
% u = gamma*w0./sqrtw2 + sqrtw2;
q = W2(d)./(W1(d).^2);
w1sqrtQ = sqrt(q).*W1(d);
v = W0(d) ./ (w1sqrtQ);
u = gamma.*v + w1sqrtQ;
end

function w = W0(d)
w = 1/2*(1+erf(d/sqrt(2)));
%w = normcdf(-d);
end

function w = W1(d)
%w = 1/2*(1+erf(d/sqrt(2)));
w = W0(d).*d + exp(-d.^2/2)/sqrt(2*pi);
end

function w = W2(d)
w = W0(d).*(1+d.^2) + d.*exp(-d.^2/2)/sqrt(2*pi);
end
