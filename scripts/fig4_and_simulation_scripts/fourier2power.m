function p = fourier2power(c)
bin = fullCube(log2(size(c,1)));
rk = sum(bin,2);
p = NaN(max(rk)+1,size(c,2));
for ii=1:size(c,2)
    p(:,ii) = accumarray(rk+1,c(:,ii).^2);
end
end

