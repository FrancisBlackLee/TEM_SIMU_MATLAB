function [anisoUs] = ExtractAtomSiteAnisoUFromCrysInfo(crysInfo)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

% find _atom_site_aniso_label
atomTypeSymbolIdx = find(strcmp('_atom_site_aniso_label', crysInfo.loopProperty(:, 1)));
atomNum = sum(~cellfun(@isempty, crysInfo.loopProperty(atomTypeSymbolIdx, :))) - 1;

% find u11 u22...
u11Idx = find(strcmp('_atom_site_aniso_U_11', crysInfo.loopProperty(:, 1)));
u22Idx = find(strcmp('_atom_site_aniso_U_22', crysInfo.loopProperty(:, 1)));
u33Idx = find(strcmp('_atom_site_aniso_U_33', crysInfo.loopProperty(:, 1)));
u23Idx = find(strcmp('_atom_site_aniso_U_23', crysInfo.loopProperty(:, 1)));
u13Idx = find(strcmp('_atom_site_aniso_U_13', crysInfo.loopProperty(:, 1)));
u12Idx = find(strcmp('_atom_site_aniso_U_12', crysInfo.loopProperty(:, 1)));

anisoUs = zeros(6, atomNum);
anisoUs(1, :) = str2double(crysInfo.loopProperty(u11Idx, 2 : atomNum + 1));
anisoUs(2, :) = str2double(crysInfo.loopProperty(u22Idx, 2 : atomNum + 1));
anisoUs(3, :) = str2double(crysInfo.loopProperty(u33Idx, 2 : atomNum + 1));
anisoUs(4, :) = str2double(crysInfo.loopProperty(u23Idx, 2 : atomNum + 1));
anisoUs(5, :) = str2double(crysInfo.loopProperty(u13Idx, 2 : atomNum + 1));
anisoUs(6, :) = str2double(crysInfo.loopProperty(u12Idx, 2 : atomNum + 1));

anisoUs = repmat(anisoUs, 1, symEquivOpNum);

end