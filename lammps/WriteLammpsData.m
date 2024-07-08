function WriteLammpsData(filename, convMat, axyz, thr)
%WriteLammpsData.m writes the structure to lammps data file.
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

if nargin == 3
    thr = 1.0e-6;
end

% if it is triclinic, use the general form
if convMat(2, 1) > thr || convMat(3, 1) > thr || convMat(3, 2) > thr
    error("The input conversion matrix is not the general form of a triclinic system");
end

fid = fopen(filename, "w");
if fid == -1
    error("Failed to open file");
else
    lx = convMat(1, 1);
    xy = convMat(1, 2);
    xz = convMat(1, 3);
    ly = convMat(2, 2);
    yz = convMat(2, 3);
    lz = convMat(3, 3);

    aTypes = axyz(1, :);
    nAtom = length(aTypes);

    uniqueATypes = uniquetol(aTypes);
    nAType = length(uniqueATypes);

    [~, sortInds] = sort(aTypes);
    axyz = axyz(:, sortInds);

    fprintf(fid, "#\n\n");
    fprintf(fid, "%d atoms\n", nAtom);
    fprintf(fid, "%d atom types\n\n", nAType);
    fprintf(fid, "0.0 %.8f xlo xhi\n", lx);
    fprintf(fid, "0.0 %.8f ylo yhi\n", ly);
    fprintf(fid, "0.0 %.8f zlo zhi\n\n", lz);
    fprintf(fid, "%.8f %.8f %.8f xy xz yz\n\n", xy, xz, yz);
    fprintf(fid, "Atom Type Labels\n\n");
    for iAType = 1 : nAType
        fprintf(fid, "%d %s\n", iAType, AtomTypeIdxToStr(uniqueATypes(iAType)));
    end
    fprintf(fid, "\n\nAtoms\n\n");
    for iAtom = 1 : nAtom
        fprintf(fid, "%d %s %.8f %.8f %.8f\n", iAtom, AtomTypeIdxToStr(axyz(1, iAtom)), ...
            axyz(2, iAtom), axyz(3, iAtom), axyz(4, iAtom));
    end

    fclose(fid);
end

end