function [a1p, a2p, a3p, b1p, b2p, b3p] = MsbtBasesTransform(a1, a2, a3, ...
    cellLengths, cellAngles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = cellLengths(1);
b = cellLengths(2);
c = cellLengths(3);

alpha = cellAngles(1);
beta = cellAngles(2);

[b1, b2, b3] = DirectToReciprocal(a1, a2, a3);
a1p = a1;
a2p = a2;
a3p = -(c / a * cosd(beta)) * a1 - (c / b * cosd(alpha)) * a2 + a3;

b1p = b1 + (c / a * cosd(beta)) * b3;
b2p = b2 + (c / b * cosd(alpha)) * b3;
b3p = b3;

end