function c = landscape2fourier(f)
% Input: a vertical vector of fitness values for a combinatorially complete set
% OR a matrix, whose columns are such vectors
% Output: the same-length vector of Fourier coefficients c
D = log2(size(f,1));
H = 1;
for i=1:D
    H = [H H;
        -H H];
end
c = H*f/(size(f,1));
end


%     V = [V/2  0*V;
%          0*V   -V];
%     H = [H H;
%          H -H];
