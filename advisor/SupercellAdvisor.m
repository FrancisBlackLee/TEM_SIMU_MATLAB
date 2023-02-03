function [superFracCoords] = SupercellAdvisor(fracCoords, zoneAxes, bases, convMat)
%SupercellAdvisor.m constructs the supercell with respect to the input unit
%cell and other related arguments.
% Input:
%   fracCoords -- fractional coordinates of the unit cell, syntax: [T; P;
%       fracX; fracY; fracZ];
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
% Output:
%   superFracCoords -- fractional coordinates of the super cell, syntax:
%       [T; P; fracX; fracY; fracZ];

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

uList = [zoneAxes.a(1), zoneAxes.b(1), zoneAxes.c(1)];
uMin = min(uList);
uMax = max(uList);
uNum = uMax - uMin + 1;

vList = [zoneAxes.a(2), zoneAxes.b(2), zoneAxes.c(2)];
vMin = min(vList);
vMax = max(vList);
vNum = vMax - vMin + 1;

wList = [zoneAxes.a(3), zoneAxes.b(3), zoneAxes.c(3)];
wMin = min(wList);
wMax = max(wList);
wNum = wMax - wMin + 1;

unitAtomNum = size(fracCoords, 2);
superAtomNum = uNum * vNum * wNum * unitAtomNum;
superFracCoords = zeros(5, superAtomNum);

for uIdx = 0 : uNum - 1
    uInc = uMin + uIdx;
    for vIdx = 0 : vNum - 1
        vInc = vMin + vIdx;
        for wIdx = 0 : wNum - 1
            wInc = wMin + wIdx;
            headIdx = uIdx * (vNum * wNum * unitAtomNum) +...
                vIdx * (wNum * unitAtomNum) + wIdx * unitAtomNum + 1;
            rearIdx = headIdx + unitAtomNum - 1;
            superFracCoords(:, headIdx : rearIdx) = fracCoords;
            superFracCoords(3, headIdx : rearIdx) = superFracCoords(3, headIdx : rearIdx) + uInc;
            superFracCoords(4, headIdx : rearIdx) = superFracCoords(4, headIdx : rearIdx) + vInc;
            superFracCoords(5, headIdx : rearIdx) = superFracCoords(5, headIdx : rearIdx) + wInc;
        end
    end
end

superCartCoords = convMat * superFracCoords(3 : 5, :);
bases.a = reshape(bases.a, 1, []);
bases.b = reshape(bases.b, 1, []);
bases.c = reshape(bases.c, 1, []);
projections = [bases.a * superCartCoords;...
               bases.b * superCartCoords;...
               bases.c * superCartCoords];

projCoeffMat = [dot(bases.a, bases.a), dot(bases.a, bases.b), dot(bases.a, bases.c);...
                dot(bases.a, bases.b), dot(bases.b, bases.b), dot(bases.b, bases.c);...
                dot(bases.a, bases.c), dot(bases.b, bases.c), dot(bases.c, bases.c)];

tmpSuperFracCoords = superFracCoords;
tmpSuperFracCoords(3 : 5, :) = projCoeffMat \ projections;
tmpSuperFracCoords(3 : 5, :) = mod(tmpSuperFracCoords(3 : 5, :), 1);
tolerance = 1.0e-6;
superFracCoords = (uniquetol(tmpSuperFracCoords', tolerance, 'ByRows', true))';

end


