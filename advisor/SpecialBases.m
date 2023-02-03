function [bases, convMat, cutoff] = SpecialBases(initConvMat, zoneAxes)
%SpecialBases.m constructs the reoriented bases of the super cell and the
%reoriented conversion matrix of the unit cell.
% Input:
%   initConvMat -- conversion matrix of the unit cell;
%   zoneAxes.a -- indices of zone axis parallel to the reoriented basis A;
%   zoneAxes.b -- indices of zone axis parallel to the reoriented basis B;
%   zoneAxes.c -- indices of zone axis parallel to the reoriented basis C,
%       same as uvw;
% Output:
%   bases.a -- reoriented lattice vector A in cartesian coordinates;
%   bases.b -- reoriented lattice vector B in cartesian coordinates;
%   bases.c -- reoriented lattice vector C in cartesian coordinates;
%   convMat -- reoriented conversion matrix with respect to the original
%       unit cell, with the basis C parallel to the view direction and the
%       projection of the basis on XOY plane parallel to X axis; it can be
%       multiplied on the left with the fractional coordinates of the
%       original unit cell to gain the cartesian coordinates;
%   cutoff -- whether the tiled super cell should be cutoff near the edges,
%       in the case, cutoff is false;

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

bases.a = CrystalIndicesToBasis(initConvMat, zoneAxes.a);
bases.b = CrystalIndicesToBasis(initConvMat, zoneAxes.b);
bases.c = CrystalIndicesToBasis(initConvMat, zoneAxes.c);

directionZ = bases.c;
directionX = bases.a;

rotMat = RotationOperator(directionZ, directionX);

convMat = rotMat * initConvMat;
bases.a = rotMat * bases.a;
bases.b = rotMat * bases.b;
bases.c = rotMat * bases.c;

cutoff = false;

end



