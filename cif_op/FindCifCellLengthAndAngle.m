function [cellLength, cellAngle] = FindCifCellLengthAndAngle(filename)
%FindCifCellLengthAndAngle() finds cell lengths and angles from the CIF
%file directly.
% Input:
%   filename -- CIF filename;
% Output:
%   cellLengths -- [a, b, c];
%   cellAngles -- [alpha, beta, gamma];

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

cellLength = zeros(1, 3);
cellAngle = zeros(1, 3);
editNum = 0;
fileID = fopen(filename, 'r');
str = ReadStrFromCif(fileID);
while ~isempty(str)
    switch str
        case '_cell_length_a'
            str = ReadStrFromCif(fileID);
            cellLength(1) = str2double(str);
            editNum = editNum + 1;
        case '_cell_length_b'
            str = ReadStrFromCif(fileID);
            cellLength(2) = str2double(str);
            editNum = editNum + 1;
        case '_cell_length_c'
            str = ReadStrFromCif(fileID);
            cellLength(3) = str2double(str);
            editNum = editNum + 1;
        case '_cell_angle_alpha'
            str = ReadStrFromCif(fileID);
            cellAngle(1) = str2double(str);
            editNum = editNum + 1;
        case '_cell_angle_beta'
            str = ReadStrFromCif(fileID);
            cellAngle(2) = str2double(str);
            editNum = editNum + 1;
        case '_cell_angle_gamma'
            str = ReadStrFromCif(fileID);
            cellAngle(3) = str2double(str);
            editNum = editNum + 1;
        otherwise
            % do nothing
    end
    
    if editNum == 6
        break;
    else
        str = ReadStrFromCif(fileID);
    end
end

fclose(fileID);

end

