classdef cflCRM_crossFeeding < communityFunctionLandscape
    properties (SetAccess='protected')
        L     % # of resources
        K     % environment (carrying capacities)
        D     % resource conversion matrix (size L by L)
        S     % # of species in pool 
        C     % resource consumption matrix; S rows
        m     % maintenance cost; length S
                
        equilibrationTime   % time to allow for equilibration
    end
    properties (Access='protected')
        propertyOfInterest % a function handle set by contructor
    end
    methods
        % constructor
        function obj = cflCRM_crossFeeding(params)
            if nargin<1
                error('Can''t construct a class instance with no params supplied.')
            end
            % params has fields: 
            %   C, D, K, m
            S = size(params.C,1);
            L = length(params.K);

            % call parent constructor
            obj@communityFunctionLandscape(S);

            switch params.propertyOfInterest
                case 'biomass'
                    obj.propertyOfInterest = @(abd,r)sum(abd);
                case 'lastResource'
                    obj.propertyOfInterest = @(abd,r)r(end);
                case 'lastSpecies'
                    obj.propertyOfInterest = @(abd,r)abd(end);
                otherwise
                    error('Property not recognized or not implemented.')
            end

            % initialize 
            obj.L = L;
            obj.S = S;
            obj.K = params.K;
            obj.D = params.D;
            obj.C = params.C;
            obj.m = params.m;
            obj.equilibrationTime = params.equilibrationTime;
        end
        function [abd, rEq] = getEquilibriumAbundance(obj, C1)
            % Returns the abundances of all S+1 species at equilibrium for the
            % specified community. Whenever C1(i)==false, necessarily abd(i)=0
            % C1 is a vector of size S. The pathogen is ALWAYS added.
            species = logical(C1);
            abd = zeros(obj.S, 1);
            [abd(species), rEq] = equilibrateEcology_CrossFeeding(obj, species);
        end
        function [abdEq, rEq] = equilibrateEcology_CrossFeeding(obj, sel)
            CC = obj.C(sel,:);
            DD = obj.D(sel,:);
            mm = obj.m(sel);
            SS = length(mm);
            r0 = obj.K(:);
            n0 = 0.1*ones(SS,1);
            x0 = [r0; n0];
            xdot = @(~,x) obj.dxdt(x,obj.K,CC,DD,mm);
            sol = ode15s(xdot, [0 obj.equilibrationTime], x0);
            xEq = sol.y(:,end);
            rEq = xEq(1:obj.L);
            abdEq = xEq(obj.L+1:end);
        end
    end
    methods (Access='protected')
        function f = getF_internal(obj,C1) % C1 = a single community signature, returns a number
            [abd, r] = getEquilibriumAbundance(obj, C1);
            f = obj.propertyOfInterest(abd,r); % Computes the property of interest
        end
    end
    methods (Static)
        function xdot = dxdt(x, K,C,D,m)
            % resources = rows
            % species = columns
            L = length(K);
            r = x(1:L)';
            n = x(L+1:end);

            r(r<0)=0; 
            n(n<0)=0;
            resourceUptake = C*(r');
            ndot = n.*( resourceUptake - m);
            rdot = K-r - r.*(n'*C) + (resourceUptake.*n)'*D;
            xdot = [rdot(:); ndot(:)];
        end

    end
end