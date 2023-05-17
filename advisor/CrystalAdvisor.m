function [superFracCoords, superConvMat, cutoff] = CrystalAdvisor(cifFilename, uvw, reduce, maxReduceNum)
%CrystalAdvisor.m generates the structure of the supercell for a crystal
%unit cell with respect to the given orientation.
% Input:
%   cifFilename -- CIF filename;
%   uvw -- unit cell orientation (view direction);
% Output:
%   superFracCoords -- fractional coordinates of the super cell, syntax:
%       [T; P; fracX; fracY; fracZ];
%   superConvMat -- conversion matrix of the supercell;
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

if nargin == 1
    uvw = [0, 0, 1];
    reduce = false;
elseif nargin == 2
    reduce = false;
elseif nargin == 3
    if reduce
        maxReduceNum = 8;
    end
elseif nargin > 4
    error("Incorrect number of input arguments");
end


crysInfo = LoadCif(cifFilename);
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
unitFracCoords = ExtractAtomSiteFromCrysInfo(crysInfo);
fullUnitFracCoords = AddEquivAtomSites(unitFracCoords);

[zoneAxes, bases, convMat, cutoff] = BasesAdvisor(cellLengths, cellAngles, uvw);
superFracCoords = SupercellAdvisor(fullUnitFracCoords, zoneAxes, bases, convMat);
superFracCoords = RemoveSymmetricAtoms(superFracCoords);
superConvMat = BasesToConvMat(bases);

if reduce
    [superFracCoords, superConvMat] = ReduceUnitCell(superFracCoords, superConvMat, maxReduceNum);
end

end


