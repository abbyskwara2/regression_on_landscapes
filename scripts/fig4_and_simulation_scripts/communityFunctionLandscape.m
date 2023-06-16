classdef (Abstract) communityFunctionLandscape < handle
    properties (SetAccess='protected', GetAccess='public')
        DIM     % dimension (number of species/mutations)
        cube  % storing the computed values of F. Sparse array with space for 2^12 entries
    end
    methods
        function obj = communityFunctionLandscape(DIM)
            if nargin<1
                DIM=10;
            end
            obj.DIM = DIM;
            if DIM>12 
                maxSize = 2^12;
            else
                maxSize = 2^DIM;
            end
            obj.cube = sparse([],[],[],2^DIM,1,maxSize);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = computeAll(obj,STOP_ON_FAILURE)
            assert(obj.DIM<=12,'computeAll is meant to be used with DIM <= 12 only.')
            if nargin<2
                STOP_ON_FAILURE = false;
            end
            C = dec2bin(0:(2^obj.DIM-1))=='1'; % all binary combinations
            f = getF(obj, C, STOP_ON_FAILURE);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = getF(obj,C,STOP_ON_FAILURE) % C = Community signature, binary vector
            if nargin<3
                STOP_ON_FAILURE = false;
            end
            assert(size(C,2)==obj.DIM,'Dimension mismatch')
            f = queryCube(obj, C);
            % find values not yet computed
            needToCompute = isnan(f);
            Csubset = C(needToCompute,:);
            NN = size(Csubset,1);
            valuesToCompute = f(needToCompute); % all NaNs
            parfor i=1:NN
            %for i=1:NN
                valuesToCompute(i) = getF_internal(obj,Csubset(i,:));
                if STOP_ON_FAILURE && isnan(valuesToCompute(i))
                    error('getF:fail','Stopping on failure as requested');
                end
            end
            storeInCube(obj, Csubset, valuesToCompute);
            f(needToCompute) = valuesToCompute;
        end
        
        % Helper functions:
        function f = queryCube(obj, C1)
            % C1 -> linear index
            ind = C2ind(obj, C1);
            % for a sparse array, absent values are 0
            % for us, 0 could actually be the value stored
            % whereas absent value means NaN
            f = full(obj.cube(ind));
            f = obj.swapNAN(f);
        end
        function ind = C2ind(obj, C)
            pow2 = 2.^(0:obj.DIM-1);
            ind = (double(C))*pow2(:)+1;
        end
        function storeInCube(obj, C, f)
            f = obj.swapNAN(f);
            ind = C2ind(obj, C);
            obj.cube(ind) = f;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function C = randomSet(obj,N,p)
            C = rand(N,obj.DIM)<p; %Random matrix with communities of approximately right size
            before = size(C,1);
            % remove duplicates
            C = unique(C,'rows'); %Removes duplicate communities
            % remove empty set 
            notNull = any(C,2); %Finds rows that aren't all zeros
            C = C(notNull, :); %Removes empty rows
            after = size(C,1);
            if after<before
                warning('Removed %d duplicate entries from regression training set', before-after);
            end
        end
    end
    methods (Static)
        function f = swapNAN(f)
            zers = f==0;
            nans = isnan(f);
            f(zers) = NaN;
            f(nans) = 0;
        end
    end
    
    methods (Abstract, Access='protected')
        getF_internal(obj,C1) % C1 = a single community signature, returns a number        
    end
end