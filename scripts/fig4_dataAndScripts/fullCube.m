function cube = fullCube(D)
assert(D<=12,'Accepting D up to 12.')
cube = dec2bin(0:(2^D-1))=='1';
end
