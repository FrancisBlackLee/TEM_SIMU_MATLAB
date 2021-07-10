function [atomSiteMat] = ExtractAtomSiteFromCrysInfo(crysInfo)
%ExtractAtomSiteFromCrysInfo() extracts atom sites (after symmetry
%operation) from the crysInfo struct output by LoadCif().
% Input:
%   crysInfo -- crysInfo struct output by LoadCif();
% Output:
%   atomSiteMat -- [atomType; occupancy; fractX; fractY; fractZ];
% NOTE:
%   fractional coords, fractX, fractY, fractZ, are not cartesian;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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

% find symmetry operations
symEquivOpIdx = find(strcmp('_symmetry_equiv_pos_as_xyz', crysInfo.loopProperty(:, 1)));
symEquivOpNum = sum(~cellfun(@isempty, crysInfo.loopProperty(symEquivOpIdx, :))) - 1;
symEquivOpMat = zeros(3, 3, symEquivOpNum);
glideMat = zeros(3, symEquivOpNum);
for i = 1 : symEquivOpNum
    [symEquivOpMat(:, :, i), glideMat(:, i)] =...
        SymmetryOperatorToMatrix(crysInfo.loopProperty{symEquivOpIdx, i + 1});
end

% find atomType
atomTypeSymbolIdx = find(strcmp('_atom_site_type_symbol', crysInfo.loopProperty(:, 1)));
atomNum = sum(~cellfun(@isempty, crysInfo.loopProperty(atomTypeSymbolIdx, :))) - 1;

% find atom occupancy
atomOccuIdx = find(strcmp('_atom_site_occupancy', crysInfo.loopProperty(:, 1)));

% find fract x, y, z
fractXIdx = find(strcmp('_atom_site_fract_x', crysInfo.loopProperty(:, 1)));
fractYIdx = find(strcmp('_atom_site_fract_y', crysInfo.loopProperty(:, 1)));
fractZIdx = find(strcmp('_atom_site_fract_z', crysInfo.loopProperty(:, 1)));

typeList = zeros(1, atomNum);
for atomIdx = 1 : atomNum
    typeList(atomIdx) = AtomTypeStrToIdx(crysInfo.loopProperty{atomTypeSymbolIdx, atomIdx + 1});
end

occuList = str2double(crysInfo.loopProperty(atomOccuIdx, 2 : atomNum + 1));

fractCoordMat = zeros(3, atomNum);
fractCoordMat(1, :) = str2double(crysInfo.loopProperty(fractXIdx, 2 : atomNum + 1));
fractCoordMat(2, :) = str2double(crysInfo.loopProperty(fractYIdx, 2 : atomNum + 1));
fractCoordMat(3, :) = str2double(crysInfo.loopProperty(fractZIdx, 2 : atomNum + 1));

atomSiteMat = [typeList; occuList; fractCoordMat];
atomSiteMat = repmat(atomSiteMat, 1, symEquivOpNum);
for i = 1 : symEquivOpNum
    head = (i - 1) * atomNum + 1;
    rear = i * atomNum;
    atomSiteMat(3 : 5, head : rear) = symEquivOpMat(:, :, i) * atomSiteMat(3 : 5, head : rear) +...
        glideMat(:, i);
end

tolerance = 1e-5;

index = find(atomSiteMat(3, :) > 1 + tolerance);
atomSiteMat(3, index) = mod(atomSiteMat(3, index), 1);
index = find(atomSiteMat(3, :) < 0 - tolerance);
atomSiteMat(3, index) = atomSiteMat(3, index) - floor(atomSiteMat(3, index));

index = find(atomSiteMat(4, :) > 1 + tolerance);
atomSiteMat(4, index) = mod(atomSiteMat(4, index), 1);
index = find(atomSiteMat(4, :) < 0 - tolerance);
atomSiteMat(4, index) = atomSiteMat(4, index) - floor(atomSiteMat(4, index));

index = find(atomSiteMat(5, :) > 1 + tolerance);
atomSiteMat(5, index) = mod(atomSiteMat(5, index), 1);
index = find(atomSiteMat(5, :) < 0 - tolerance);
atomSiteMat(5, index) = atomSiteMat(5, index) - floor(atomSiteMat(5, index));

atomSiteMat = atomSiteMat';
atomSiteMat = uniquetol(atomSiteMat, tolerance, 'ByRows', true);
atomSiteMat = atomSiteMat';

end

