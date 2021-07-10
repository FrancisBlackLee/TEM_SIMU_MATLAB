function [typeList, xyzList] = ReadCrystalMakerXyz(filename)
%ReadCrystalMakerXyz() reads crystal coordinate data from xyz file export
%from CrystalMaker.
% Input:
%   filename -- filename of the xyz file;
% Output:
%   typeList -- list of the atomic types, 1 by atomNum array;
%   xyzList -- list of the atomic xyz coordinates, 3 by atomNum matrix;

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

filePtr = fopen(filename, 'r');
if filePtr == -1
    typeList = [];
    xyzList = [];
    msgbox('Failed to open file!');
else
    atomNum = fscanf(filePtr, '%d', 1);
    typeList = zeros(1, atomNum);
    xyzList = zeros(3, atomNum);
    for atomIdx = 1 : atomNum
        tmpTypeStr = fscanf(filePtr, '%s', 1);
        typeList(atomIdx) = AtomTypeStrToIdx(tmpTypeStr);
        xyzList(1, atomIdx) = fscanf(filePtr, '%f', 1);
        xyzList(2, atomIdx) = fscanf(filePtr, '%f', 1);
        xyzList(3, atomIdx) = fscanf(filePtr, '%f', 1);
    end

    fclose(filePtr);
end

end

