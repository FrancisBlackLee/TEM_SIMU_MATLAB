function [zoneAxes, bases, convMat] = NormalBasesAdvisor(cellLengths,...
    cellAngles, uvw)
%NormalBasesAdvisor.m provides under the preset rule the advised reoriented
%zone axes, lattice vectors (bases) and the corresponding conversion matrix
%with respect to the original unit cell. This function can be applied to
%nonspecial orientations of crystal systems except the cubic system, and
%the function is suitable to all orientations of the cubic system.
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
zoneAxes.c = uvw;
bases.c = CrystalIndicesToBasis(initConvMat, zoneAxes.c);

% ReshaperUV
[zoneAxes.a, zoneAxes.b] = ReshaperUV;
bases.a = CrystalIndicesToBasis(initConvMat, zoneAxes.a);
bases.b = CrystalIndicesToBasis(initConvMat, zoneAxes.b);
volume = abs(dot(bases.a, cross(bases.b, bases.c)));

% ReshaperVW
[tmpIndicesA, tmpIndicesB] = ReshaperVW;
tmpBasisA = CrystalIndicesToBasis(initConvMat, tmpIndicesA);
tmpBasisB = CrystalIndicesToBasis(initConvMat, tmpIndicesB);
tmpVolume = abs(dot(tmpBasisA, cross(tmpBasisB, bases.c)));
if tmpVolume < volume
    zoneAxes.a = tmpIndicesA;
    zoneAxes.b = tmpIndicesB;
    bases.a = tmpBasisA;
    bases.b = tmpBasisB;
    volume = tmpVolume;
end

% ReshaperWU
[tmpIndicesA, tmpIndicesB] = ReshaperWU;
tmpBasisA = CrystalIndicesToBasis(initConvMat, tmpIndicesA);
tmpBasisB = CrystalIndicesToBasis(initConvMat, tmpIndicesB);
tmpVolume = abs(dot(tmpBasisA, cross(tmpBasisB, bases.c)));
if tmpVolume < volume
    zoneAxes.a = tmpIndicesA;
    zoneAxes.b = tmpIndicesB;
    bases.a = tmpBasisA;
    bases.b = tmpBasisB;
    volume = tmpVolume;
end

directionZ = bases.c;
directionX = bases.a - dot(bases.a, bases.c) * bases.c / norm(bases.c)^2;
rotMat = RotationOperator(directionZ, directionX);

convMat = rotMat * initConvMat;
bases.a = rotMat * bases.a;
bases.b = rotMat * bases.b;
bases.c = rotMat * bases.c;

% nested functions
    function [indicesA, indicesB] = ReshaperUV
        indicesA = [v, -u, 0];
        if ~any(indicesA)
            indicesA = [1, 0, 0];
            indicesB = [0, 1, 0];
        else
            indicesB = [u * w, v * w, -u^2 - v^2];

            % minimize indicesA by the greatest common divisor of its non-zero elements
            vecAGcd = gcd(u, v);
            if vecAGcd == 0
                indicesA = indicesA / max([u, v]);
            else
                indicesA = indicesA / vecAGcd;
            end

            % minimize indicesB by the greatest common divisor of its non-zero elements
            vecBNonZeros = indicesB(indicesB ~= 0);
            indicesB = indicesB / double(gcd(sym(vecBNonZeros)));
        end
    end

    function [indicesA, indicesB] = ReshaperVW
        indicesA = [0, w, -v];
        if ~any(indicesA)
            indicesA = [0, 1, 0];
            indicesB = [0, 0, 1];
        else
            indicesB = [-v^2 - w^2, u * v, u * w];
            
            % minimize indicesA by the greatest common divisor of its non-zero elements
            vecAGcd = gcd(v, w);
            if vecAGcd == 0
                indicesA = indicesA / max([v, w]);
            else
                indicesA = indicesA / vecAGcd;
            end
            
            % minimize indicesB by the greatest common divisor of its non-zero elements
            vecBNonZeros = indicesB(indicesB ~= 0);
            indicesB = indicesB / double(gcd(sym(vecBNonZeros)));
        end
    end

    function [indicesA, indicesB] = ReshaperWU
        indicesA = [w, 0, -u];
        if ~any(indicesA)
            indicesA = [1, 0, 0];
            indicesB = [0, 0, 1];
        else
            indicesB = [-u * v, w^2 + u^2, -v * w];
            
            % minimize indicesA by the greatest common divisor of its non-zero elements
            vecAGcd = gcd(w, u);
            if vecAGcd == 0
                indicesA = indicesA / max([w, u]);
            else
                indicesA = indicesA / vecAGcd;
            end
            
            % minimize indicesB by the greatest common divisor of its non-zero elements
            vecBNonZeros = indicesB(indicesB ~= 0);
            indicesB = indicesB / double(gcd(sym(vecBNonZeros)));
        end
    end

end



