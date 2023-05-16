function [b1, b2, b3] = DirectToReciprocal(a1, a2, a3)
%DirectToReciprocal calculates the reciprocal lattice vectors using the
%direct lattice vectors.
%   a1, a2, a3 -- direct lattice vectors;
%   b1, b2, b3 -- reciprocal lattice vectors;

vc = dot(a1, cross(a2, a3));
b1 = 2 * pi * cross(a2, a3) / vc;
b2 = 2 * pi * cross(a3, a1) / vc;
b3 = 2 * pi * cross(a1, a2) / vc;

end