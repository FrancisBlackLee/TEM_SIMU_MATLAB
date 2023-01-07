function [typeList, xyzList, label] = ReadVestaXyz(filename)
%ReadCrystalMakerXyz() reads crystal coordinate data from xyz file export
%from VESTA.
% Input:
%   filename -- filename of the xyz file;
% Output:
%   typeList -- list of the atomic types, 1 by atomNum array;
%   xyzList -- list of the atomic xyz coordinates, 3 by atomNum matrix;

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

if isfile(filename)
    lines = readlines(filename);
    atomNum = str2num(lines{1});
    label = lines{2};
    typeList = zeros(1, atomNum);
    xyzList = zeros(3, atomNum);
    for atomIdx = 1 : atomNum
        names = split(lines{atomIdx + 2});
        if isempty(names{1})
            names(1) = [];
        end
        typeList(atomIdx) = AtomTypeStrToIdx(names{1});
        xyzList(1, atomIdx) = str2double(names{2});
        xyzList(2, atomIdx) = str2double(names{3});
        xyzList(3, atomIdx) = str2double(names{4});
    end
else
    error('File does not exist');
end

end



