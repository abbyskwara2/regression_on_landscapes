function f = fourier2landscape(c)
% Input: a vector of fitness values for a combinatorially complete set
% Output: the same-length vector of Fourier coefficients
D = log2(length(c));
M = 1;
for i=1:D
    M = [M -M;
         M  M];
end
f = M*c(:);
end


% for i=1:D
%     M = [M -M/2;
%          M  M/2];
% end
% f = M*c(:);
