function SimpleCifGenerator(filename, convMat, fracCoords)
%SimpleCifGenerator.m automatically generates a simple cif for given
%conversion matrix and fractional coordinates.
%   filename -- filename of cif;
%   convMat -- conversion matrix;
%   fracCoords -- fractional coordinates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2022  Francis Black Lee (Li Xian)

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

bases = ConvMatToBases(convMat);
% cellLengths: [a, b, c]
a = norm(bases.a);
b = norm(bases.b);
c = norm(bases.c);

vc = abs(dot(bases.a, cross(bases.b, bases.c)));

% cellAngles: [alpha (b^c), beta (a^c), gamma (a^b)]
alpha = acosd(dot(bases.b, bases.c) / (b * c));
beta = acosd(dot(bases.a, bases.c) / (a * c));
gamma = acosd(dot(bases.a, bases.b) / (a * b));

newFracCoords = RemoveSymmetricAtoms(fracCoords);
atomNum = size(newFracCoords, 2);
% sort by atomic number Z
[~, atomTypeOrder] = sort(newFracCoords(1, :), 'ascend');
newFracCoords = newFracCoords(:, atomTypeOrder);

% write info to dest CIF
fileID = fopen(filename, 'w');

fprintf(fileID, '# Simple CIF generated by TEM_SIMU_MATLAB\n');
fprintf(fileID, '_symmetry_space_group_name_H-M    ''P 1''\n');
fprintf(fileID, '_cell_length_a    %.6f\n', a);
fprintf(fileID, '_cell_length_b    %.6f\n', b);
fprintf(fileID, '_cell_length_c    %.6f\n', c);
fprintf(fileID, '_cell_angle_alpha    %.6f\n', alpha);
fprintf(fileID, '_cell_angle_beta    %.6f\n', beta);
fprintf(fileID, '_cell_angle_gamma    %.6f\n', gamma);
fprintf(fileID, '_symmetry_Int_Tables_number    1\n');
fprintf(fileID, '_cell_volume    %.6f\n', vc);
fprintf(fileID, 'loop_\n');
fprintf(fileID, ' _symmetry_equiv_pos_site_id\n');
fprintf(fileID, ' _symmetry_equiv_pos_as_xyz\n');
fprintf(fileID, '  1  ''x, y, z''\n');
fprintf(fileID, 'loop_\n');
fprintf(fileID, ' _atom_site_type_symbol\n');
fprintf(fileID, ' _atom_site_label\n');
fprintf(fileID, ' _atom_site_symmetry_multiplicity\n');
fprintf(fileID, ' _atom_site_fract_x\n');
fprintf(fileID, ' _atom_site_fract_y\n');
fprintf(fileID, ' _atom_site_fract_z\n');
fprintf(fileID, ' _atom_site_occupancy\n');

% write atoms:
tmpType = newFracCoords(1, 1);
tmpAtomSymbol = AtomTypeIdxToStr(tmpType);
atomCount = 0;
for atomIdx = 1 : atomNum
    if newFracCoords(1, atomIdx) ~= tmpType
        tmpType = newFracCoords(1, atomIdx);
        tmpAtomSymbol = AtomTypeIdxToStr(tmpType);
    end
    tmpAtomLabel = [tmpAtomSymbol, num2str(atomCount)];
    fprintf(fileID, '  %s  %s  1  %.6f  %.6f  %.6f  %.6f\n',...
        tmpAtomSymbol, tmpAtomLabel, newFracCoords(3, atomIdx),...
        newFracCoords(4, atomIdx), newFracCoords(5, atomIdx),...
        newFracCoords(2, atomIdx));
    atomCount = atomCount + 1;
end

fclose(fileID);

end
