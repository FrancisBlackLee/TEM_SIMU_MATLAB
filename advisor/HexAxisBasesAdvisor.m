function [zoneAxes, bases, convMat, cutoff] = HexAxisBasesAdvisor(cellLengths,...
    cellAngles, uvw)
%HexAxisBasesAdvisor.m provides under the preset rule the advised
%reoriented zone axes, lattice vectors (bases) and the corresponding
%conversion matrix of the unit cell for trigonal and hexagonal crystal
%systems in hexagonal axes.
% Input:
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases b and c);
%       2 for cell angle beta (between bases a and c) and
%       3 for cell angle gamma (between bases a and b);
%   uvw -- indices of zone axis parallel to the view (Z) direction;
% Output:
%   zoneAxes.a -- indices of zone axis parallel to the reoriented basis A;
%   zoneAxes.b -- indices of zone axis parallel to the reoriented basis B;
%   zoneAxes.c -- indices of zone axis parallel to the reoriented basis C,
%       same as uvw;
%   bases.a -- reoriented lattice vector A in cartesian coordinates;
%   bases.b -- reoriented lattice vector B in cartesian coordinates;
%   bases.c -- reoriented lattice vector C in cartesian coordinates;
%   convMat -- reoriented conversion matrix with respect to the original
%       unit cell, with the basis C parallel to the view direction and the
%       projection of the basis on XOY plane parallel to X axis; it can be
%       multiplied on the left with the fractional coordinates of the
%       original unit cell to gain the cartesian coordinates.
%   cutoff -- whether the tiled super cell should be cutoff near the edges,
%       if uvw is a special orientation, under which the periodic boundary
%       conditions is satisfied, cutoff is false; otherwise cutoff should 
%       be true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~any(uvw, 'all')
    error('Invalid input of uvw');
end

uvwNonzeros = uvw(uvw ~= 0);
uvw = uvw / double(gcd(sym(uvwNonzeros)));

u = uvw(1);
v = uvw(2);
w = uvw(3);

uvw = reshape(uvw, 1, []);
initConvMat = ConversionMatrix(cellLengths, cellAngles);

% special case 1
if isequal(abs(uvw), [0, 0, 1])
    zoneAxes.a = sign(w) * [1, -1, 0];
    zoneAxes.b = [1, 1, 0];
    zoneAxes.c = uvw;
    [bases, convMat, cutoff] = SpecialBases(initConvMat, zoneAxes);
% special case 2
elseif isequal(uvw, [1, 2, 0]) || isequal(uvw, [-1, -2, 0])
    zoneAxes.a = sign(u) * [-1, 0, 0];
    zoneAxes.b = [0, 0, 1];
    zoneAxes.c = uvw;
    [bases, convMat, cutoff] = SpecialBases(initConvMat, zoneAxes);
% special case 3
elseif isequal(uvw, [2, 1, 0]) || isequal(uvw, [-2, -1, 0])
    zoneAxes.a = sign(u) * [0, 1, 0];
    zoneAxes.b = [0, 0, 1];
    zoneAxes.c = uvw;
    [bases, convMat, cutoff] = SpecialBases(initConvMat, zoneAxes);
% special case 4
elseif (w == 0) && (u ~= 2 * v) && (v ~= 2 * u)
    ua = u - 2 * v;
    va = 2 * u - v;
    uavaGcd = gcd(ua, va);
    ua = ua / uavaGcd;
    va = va / uavaGcd;
    zoneAxes.a = [ua, va, 0];
    zoneAxes.b = [0, 0, 1];
    zoneAxes.c = uvw;
    [bases, convMat, cutoff] = SpecialBases(initConvMat, zoneAxes);
else
    % go to NormalBaseAdvisor.m
    [zoneAxes, bases, convMat] = NormalBasesAdvisor(cellLengths, cellAngles, uvw);
    cutoff = true;
end

end



