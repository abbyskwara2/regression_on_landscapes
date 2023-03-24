classdef cflLotkaVolterra < communityFunctionLandscape
    properties (SetAccess='protected')
        S % # of species in pool
        inters % Interaction matrix, S x S
        growth % Growth rates, S
        K % Carrying capacities, S
        abd0 %Starting abundances S

        Teq % Equilibration time
    end
    properties (Access='protected')
        propertyOfInterest % a function handle set by contructor
        solveropt % setting ode options object takes time, so only do this once
    end

    methods 
        % This constructor builds our community & calculates abundance of
        % the pathogen. 
        function obj = cflLotkaVolterra(params)
            if nargin<1 %Error handling if no parameters provided
                error('Can''t construct a class instance with no params supplied.')
            end
            % params contain the info needed for LV random model
            % S = number of species in the pool
            % inters = interaction matrix (S x S)
            % growth = growth vector, growth rate for each species
            % K = carrying capacity vector, carry. cap. for each species
            % abd0 = initial abundances vector, initial abd. for each
            %           species
            % Teq = equilibration time
            obj@communityFunctionLandscape(params.S);
            
            %% Setting up initial data
            if isa(params.propertyOfInterest,'function_handle')
                obj.propertyOfInterest = params.propertyOfInterest;
            else
                switch params.propertyOfInterest
                    case 'biomass'
                        obj.propertyOfInterest = @(abd)sum(abd);
                    otherwise
                        error('Property not recognized or not implemented.')
                end
            end
            %Pass params to object
            obj.S = params.S;
            obj.Teq = params.Teq;
            obj.inters = params.inters;
            obj.K = params.K(:); 
            obj.growth = params.growth(:);
            obj.abd0 = params.abd0(:);
            obj.solveropt = odeset('Events',@obj.detectDivergence);
        end

        function abd = getAbd(obj,C1) 
            % wrapper for getF_internal that can be accessed from outside
            % C1 = a single binary community signature
            if any(C1)
                [abd, failure] = evolve_lv_model(obj, C1);
            else
                abd = zeros(obj.S,1);
                failure = false;
            end
            if failure
                abd = NaN(size(abd));
            end
        end

    end
    methods (Access = 'protected')
        %This function evolves our community w/ LV dynamics 
        function [abds, failure] = evolve_lv_model(obj,C1)
            % C1 describes species other than the pathogen. 
            % The pathogen is always added.
            select = C1(:); % a column vector of species present

            % Restrict the full LV system to the species that are present
            n0 = obj.abd0(select);
            RR = obj.growth(select);
            AA = obj.inters(select, select);
            KK = obj.K(select);
            function dn = dndt(~,n)
                dn = RR./KK .* n .* (KK - n - AA*n); % Following Eq (1) in Bunin 1607.04734
                % enforce that if n<0, dn can only be positive
                dn(n<0 & dn<0) = 0;
            end

            [~, n, te, ye, ie] = ode15s(@dndt, [0 obj.Teq], n0, obj.solveropt);
            failure = ~isempty(te); 
            abds = zeros(obj.S,1);
            abds(select) = n(end,:);
        end    
        function [f, abd] = getF_internal(obj,C1) 
            % C1 = a single binary community signature, returns a number
            if any(C1)
                [abd, failure] = evolve_lv_model(obj, C1);
            else
                abd = zeros(obj.S,1);
                failure = false;
            end
            if failure
                f = NaN;
            else
                f = obj.propertyOfInterest(abd); %Computes the property of interest
            end
        end
    end
    methods (Static)
        function [value,isterminal,direction] = detectDivergence(~,y)
            value = max(y) - 1e6; % Detect if any abundance crosses threshold
            isterminal = 1;   % Stop the integration
            direction = [];   % 
        end
    end
end
