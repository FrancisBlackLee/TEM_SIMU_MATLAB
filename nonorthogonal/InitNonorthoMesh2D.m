function [xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2)
%InitNonorthoMesh2D generates the real space and reciprocal meshes provided
%the 2D real space lattice vectors, repeating numbers and sampling along
%each lattice vector.
% [xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2)
%   a1 -- lattice vector 1;
%   a2 -- lattice vector 2;
%   na1 -- repeating number along lattice vector 1;
%   na2 -- repeating number along lattice vector 2;
%   n1 -- sampling number along lattice vector 1;
%   n2 -- sampling number along lattice vector 2;
%   xMesh -- x mesh
%   yMesh -- y mesh
%   fxMesh -- fx mesh
%   fyMesh -- fy mesh

% da1 = a1 * na1 / n1;
% da2 = a2 * na2 / n2;
% da3 = [0, 0, 1];
% 
% [b1, b2, ~] = DirectToReciprocal(da1, da2, da3);

a3 = [0, 0, 1];
[b1, b2, ~] = DirectToReciprocal(a1, a2, a3);
b1 = b1 * n1 / na1 / (2 * pi);
b2 = b2 * n2 / na2 / (2 * pi);

% real space:
gridA1 = InitAxis(na1, n1);
gridA2 = InitAxis(na2, n2);

[meshA1, meshA2] = meshgrid(gridA1, gridA2);
xMesh = a1(1) * meshA1 + a2(1) * meshA2;
yMesh = a1(2) * meshA1 + a2(2) * meshA2;

% reciprocal space:
gridB1 = InitAxis(1, n1);
gridB2 = InitAxis(1, n2);

[meshB1, meshB2] = meshgrid(gridB1, gridB2);
fxMesh = b1(1) * meshB1 + b2(1) * meshB2;
fyMesh = b1(2) * meshB1 + b2(2) * meshB2;

end