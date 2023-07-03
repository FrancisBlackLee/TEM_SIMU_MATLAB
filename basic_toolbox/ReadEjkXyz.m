function [lattConsts, typeCoords, wobbles] = ReadEjkXyz(filename)
%ReadEjkXyz reads an xyz file of E. J. Kirkland's format.
% Input:
%   filename -- filename of the input ejk xyz file;
% Output:
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

fileID = fopen(filename);
if fileID == -1
    error('Failed to open file!');
else
    commentLine = fgetl(fileID);
    lattConsts(1) = fscanf(fileID, '%f', 1);
    lattConsts(2) = fscanf(fileID, '%f', 1);
    lattConsts(3) = fscanf(fileID, '%f', 1);
    tmpType = fscanf(fileID, '%d', 1);
    atomCount = 0;
    while tmpType ~= -1
        atomCount = atomCount + 1;
        typeCoords(1, atomCount) = tmpType;
        typeCoords(3, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(4, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(5, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(2, atomCount) = fscanf(fileID, '%f', 1);
        wobbles(atomCount) = fscanf(fileID, '%f', 1);
        
        tmpType = fscanf(fileID, '%d', 1);
    end
    fclose(fileID);
end

end
