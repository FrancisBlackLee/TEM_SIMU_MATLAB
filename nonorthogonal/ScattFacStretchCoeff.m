function stretchCoeff = ScattFacStretchCoeff(a1, a2, na1, na2, n1, n2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

t = [a1(1), a2(1); a1(2), a2(2)];
d1 = na1 / n1;
d2 = na2 / n2;

stretchCoeff = abs(1 / (det(t) * d1 * d2));

end