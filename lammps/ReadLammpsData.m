function [convMat, ids, axyz] = ReadLammpsData(filename)
%WriteLammpsData.m reads the structure from lammps data file.
%   filename -- filename of the lammps data file
%   convMat -- conversion matrix, refer to Lammps document for more
%       details;
%   axyz -- atomic types and coordinates, 4-by-N matrix, row 1: atomic
%       types; row 2 ~ row 4: cartesian coordinates.

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

fid = fopen(filename, "r");
if fid == -1
    error("Failed to open file");
else
    while ~feof(fid)
        tmpLine = fgetl(fid);
        
        if contains(tmpLine, 'atoms')
            nAtom = sscanf(tmpLine, '%d');
        elseif contains(tmpLine, 'atom') && contains(tmpLine, 'types')
            nType = sscanf(tmpLine, '%d');
        elseif contains(tmpLine, 'xlo') && contains(tmpLine, 'xhi')
            xlo_xhi = sscanf(tmpLine, '%f');
            lx = xlo_xhi(2) - xlo_xhi(1);
        elseif contains(tmpLine, 'ylo') && contains(tmpLine, 'yhi')
            ylo_yhi = sscanf(tmpLine, '%f');
            ly = ylo_yhi(2) - ylo_yhi(1);
        elseif contains(tmpLine, 'zlo') && contains(tmpLine, 'zhi')
            zlo_zhi = sscanf(tmpLine, '%f');
            lz = zlo_zhi(2) - zlo_zhi(1);
        elseif contains(tmpLine, 'xy') && contains(tmpLine, 'xz') && contains(tmpLine, 'yz')
            xy_xz_yz = sscanf(tmpLine, '%f');
            xy = xy_xz_yz(1);
            xz = xy_xz_yz(2);
            yz = xy_xz_yz(3);
        elseif contains(tmpLine, 'Atoms')
            break;
        end
    end

    convMat = zeros(3);
    convMat(1, 1) = lx;
    convMat(1, 2) = xy;
    convMat(1, 3) = xz;
    convMat(2, 2) = ly;
    convMat(2, 3) = yz;
    convMat(3, 3) = lz;

    ids = zeros(1, nAtom);
    axyz = zeros(4, nAtom);
    atomCount = 0;

    while ~feof(fid)
        tmpLine = fgetl(fid);

        if ~isempty(tmpLine)
            atomCount = atomCount + 1;
            strs = split(tmpLine);
            ids(atomCount) = str2double(strs{1});
            axyz(1, atomCount) = AtomTypeStrToIdx(strs{2});
            axyz(2, atomCount) = str2double(strs{3});
            axyz(3, atomCount) = str2double(strs{4});
            axyz(4, atomCount) = str2double(strs{5});
        end
    end

    fclose(fid);
end

end