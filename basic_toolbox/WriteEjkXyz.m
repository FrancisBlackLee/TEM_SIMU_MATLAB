function WriteEjkXyz(filename, lattConsts, typeCoords, wobbles)
%WriteEjkXyz writes an xyz file of E. J. Kirkland's format.
% Input:
%   filename -- filename of the input ejk xyz file;
%   lattConsts -- lattice constants in angstrom;
%   typeCoords -- a matrix describing the unit cell, the frist row denotes
%       the atomic number, second row the atomic occupancy, third to fifth
%       rows the atomic cartesian coordinates;
%   wobbles -- mean square root of atomic displacements.

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

fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file!');
else
    fprintf(fid, 'Created by TEM_SIMU_MATLAB\n');
    fprintf(fid, '\t%.6f\t%.6f\t%.6f\n', lattConsts(1), lattConsts(2),...
        lattConsts(3));
    atomNum = size(typeCoords, 2);
    for atomIdx = 1 : atomNum
        fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',...
            typeCoords(1, atomIdx), typeCoords(3, atomIdx), typeCoords(4, atomIdx),...
            typeCoords(5, atomIdx), typeCoords(2, atomIdx), wobbles(atomIdx));
    end
    fprintf(fid, '-1\n');

    fclose(fid);
end

end