function WriteVestaXyz(filename, typeList, xyzList, label)
%WriteCrystalMakerXyz() writes crystal coordinate data as VESTA xyz file.
% Input:
%   filename -- filename of the xyz file;
%   typeList -- list of the atomic types, 1 by atomNum array;
%   xyzList -- list of the atomic xyz coordinates, 3 by atomNum matrix;

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

filePtr = fopen(filename, 'w');
if filePtr == -1
    error('Failed to open file!');
else
    atomNum = length(typeList);
    fprintf(filePtr, '%d\n', atomNum);
    fprintf(filePtr, '%s\n', label);
    for atomIdx = 1 : atomNum
        fprintf(filePtr, '%s  %f  %f  %f\n',...
            AtomTypeIdxToStr(typeList(atomIdx)),...
            xyzList(1, atomIdx),...
            xyzList(2, atomIdx),...
            xyzList(3, atomIdx));
    end

    fclose(filePtr);
end

end


