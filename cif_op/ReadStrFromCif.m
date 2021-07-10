function [str] = ReadStrFromCif(fileID)
%ReadStrFromCif() read a valid string from the CIF file.
% Input:
%   fileID -- file identity/handle/pointer of the file;
% Output:
%   str -- string;

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

str = fscanf(fileID, '%s', 1);

if ~isempty(str)
    if str(1) == ''''
        if str(end) ~= ''''
            tmpStr = fscanf(fileID, '%s', 1);
            while tmpStr(end) ~= ''''
                str = [str, tmpStr];
                tmpStr = fscanf(fileID, '%s', 1);
            end
            str = [str, tmpStr];
        end
        str(1) = [];
        str(end) = [];
    end
    
    % check if parentheses exist, if they do, delete them and contents
    % inside
    if contains(str, '(') && contains(str, ')')
        leftParenIdx = strfind(str, '(');
        rightParenIdx = strfind(str, ')');
        str(leftParenIdx : rightParenIdx) = [];
    end
end

end

