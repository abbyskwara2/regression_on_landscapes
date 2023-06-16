function [circuit, G, iter] = crossFeedingCircuitInstance(S,L,lambda,seed)
% Circuit = what you eat (C) and what you secrete (D)
% Define species with 1 resource in, 1 resource out (or possibly none out)
% Resource in = purely random
% Resource out = prob. lambda to be resource i+1, otherwise all equiprob
% 	L resources in cascade
% 	S species, S>=L
% A circuit is valid if last resource can be produced from the first one

assert(S>=L);

if S == L
    % the pure linear pathway is the only choice
    circuit.C = eye(S); 
    circuit.D = diag(ones(S-1,1),1);
    circuit.S = S;
    circuit.L = L;
    iter = 0;
else
    rng(seed);
    keepLooking = true;
    iter = 0;
    while keepLooking 
        circuit = getCircuitInstance(S,L,lambda);
        %imagesc(circuit.C-circuit.D);

        [valid, G] = isValidCircuit(circuit);
        keepLooking = ~valid;
        iter = iter+1;
    end
end
end


function circuit = getCircuitInstance(S,L,lambda)
% 	Consume & secrete: Two samples from a Gaussian
% 	Computed as i0+-di
%   i0 distributed uniformly between rho/2 and L-rho/2
%   di is from a Gaussian of width rho/2.
from = randi(L, [S,1]);
% pick TO by picking number of steps
prob = [lambda, (1-lambda)*ones(1,L-1)/(L-1)];
steps = randsample(L,S,true,prob);
to = mod(from+steps, L+1);
% to=0 means "secrete nothing"

[from, ind] = sort(from);
to = to(ind);

circuit.from = from;
circuit.to = to;
circuit.to(to==0) = L+1;
circuit.G = digraph(circuit.from, circuit.to, [], L+1);

idx = (1:S)';
C = zeros(S,L);
D = zeros(S,L);
C(sub2ind(size(C),idx,from))=1;

idx = idx(to>0);
to = to(to>0);
D(sub2ind(size(D),idx,to))=1;

circuit.C = C;
circuit.D = D;
circuit.S = S;
circuit.L = L;
end


% a circuit is "valid" if last resource is reachable from the first
function [valid, G] = isValidCircuit(circuit)
    adjacencyMatrix = false(circuit.L, circuit.L);
    for i=1:size(circuit.C,1)
        from = circuit.C(i,:)>0;
        to = circuit.D(i,:)>0;
        adjacencyMatrix(from,to)=true;
    end
    G = digraph(adjacencyMatrix);
    path = shortestpath(G,1,circuit.L);
    valid = ~isempty(path);
end