function [atomSites] = ExtractAtomSiteFromCrysInfo(crysInfo)
%ExtractAtomSiteFromCrysInfo() extracts atom sites (after symmetry
%operation) from the crysInfo struct output by LoadCif().
% Input:
%   crysInfo -- crysInfo struct output by LoadCif();
% Output:
%   atomSiteMat -- [atomType; occupancy; fractX; fractY; fractZ];
% NOTE:
%   fractional coords, fractX, fractY, fractZ, are not cartesian;

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

% find symmetry operations
symEquivOpIdx = find(strcmp('_symmetry_equiv_pos_as_xyz', crysInfo.loopProperty(:, 1)));
% find the alternative
if isempty(symEquivOpIdx)
    symEquivOpIdx = find(strcmp('_space_group_symop_operation_xyz', crysInfo.loopProperty(:, 1)));
end

if isempty(symEquivOpIdx)
    loopPropertyNum = size(crysInfo.loopProperty, 1);
    symEquivOpIdx = loopPropertyNum + 1;
    crysInfo.loopProperty{symEquivOpIdx, 1} = '_symmetry_equiv_pos_as_xyz';
    crysInfo.loopProperty{symEquivOpIdx, 2} = 'x,y,z';
end
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

atomSites = [typeList; occuList; fractCoordMat];
atomSites = repmat(atomSites, 1, symEquivOpNum);
for i = 1 : symEquivOpNum
    head = (i - 1) * atomNum + 1;
    rear = i * atomNum;
    atomSites(3 : 5, head : rear) = symEquivOpMat(:, :, i) * atomSites(3 : 5, head : rear) +...
        glideMat(:, i);
end

tolerance = 1e-5;

index = find(atomSites(3, :) > 1 + tolerance);
atomSites(3, index) = mod(atomSites(3, index), 1);
index = find(atomSites(3, :) < 0 - tolerance);
atomSites(3, index) = atomSites(3, index) - floor(atomSites(3, index));

index = find(atomSites(4, :) > 1 + tolerance);
atomSites(4, index) = mod(atomSites(4, index), 1);
index = find(atomSites(4, :) < 0 - tolerance);
atomSites(4, index) = atomSites(4, index) - floor(atomSites(4, index));

index = find(atomSites(5, :) > 1 + tolerance);
atomSites(5, index) = mod(atomSites(5, index), 1);
index = find(atomSites(5, :) < 0 - tolerance);
atomSites(5, index) = atomSites(5, index) - floor(atomSites(5, index));

atomSites = atomSites';
atomSites = uniquetol(atomSites, tolerance, 'ByRows', true);
atomSites = atomSites';

end



