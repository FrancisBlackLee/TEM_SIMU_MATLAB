function [strCellArray, lineType, canDelete] = CifLineParser(textLine, lastLoopObj)
%CifLineParser() is the parser for the text obtained using CIF line scan.
% Input:
%   textLine -- CIF text line;
% Output:
%   strCellArray -- string cell array obtained from the text line;
%   lineType -- type of the text line classified inside the function;
%   canDelete - whether this line can be deleted, e.g., a comment line, a
%       loop commad, 'true' for 'can delete', 'false' for 'cannot delete';

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

tightTextLine = TightenCifTextLine(textLine);
% disp(tightTextLine);
strCellArray = strsplit(tightTextLine);

strNum = numel(strCellArray);
if strNum == 1
    if ~isempty(strCellArray{1})
        if strCellArray{1}(1) == '#'
            lineType = 'Comment';
            canDelete = 'true';
        elseif strCellArray{1}(1) == '_'
            lineType = 'LoopProperty';
            canDelete = 'false';
        elseif strcmp(strCellArray{1}, 'loop_')
            lineType = 'LoopCommand';
            canDelete = 'true';
        else
            if ~strcmp(lastLoopObj, 'None')
                lineType = 'LoopValues';
                canDelete = 'false';
            else
                lineType = 'Undef';
                canDelete = 'true';
            end
        end
    else
        lineType = 'WhiteLine';
        canDelete = 'true';
    end
elseif strNum == 2
    if strCellArray{1}(1) == '#'
        lineType = 'Comment';
        canDelete = 'true';
    elseif strCellArray{1}(1) == '_'
        lineType = 'ValuedProperty';
        canDelete = 'false';
    else
        lineType = 'LoopValues';
        canDelete = 'false';
    end
elseif strNum > 2
    if strCellArray{1}(1) == '#'
        lineType = 'Comment';
        canDelete = 'true';
    else
        lineType = 'LoopValues';
        canDelete = 'false';
    end
else
    lineType = 'Undef';
    canDelete = 'true';
end

end

