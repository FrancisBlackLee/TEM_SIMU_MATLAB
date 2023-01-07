function [crysInfo] = LoadCif(filename)
% LoadCif() loads a cif file and reads the basic information from the cif
% file for other functions of TEM_SIMU_MATLAB.

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

% NOTE:
%   str starts with 'data' is the simple description of the cif file;
%   str starts with underscore '_' is a property;
%   str 'loop_' means loop over the properties under it;

% first load all the text from the destination cif file

fileID = fopen(filename, 'r');

valuedPropertyNum = 0;
loopObj = 'None';
loopPropertyHead = 0;
loopPropertyNum = 0;
loopValueColNum = 2;

% disp('CIF content');
tmpTextLine = ReadTextLine;
while ischar(tmpTextLine)
    [strCellArray, lineType, canDelete] = CifLineParser(tmpTextLine, loopObj);
%     disp(strCellArray);
    if strcmp(canDelete, 'false')
        if strcmp(lineType, 'ValuedProperty')
            valuedPropertyNum = valuedPropertyNum + 1;
            crysInfo.valuedProperty{valuedPropertyNum, 1} = strCellArray{1};
            crysInfo.valuedProperty{valuedPropertyNum, 2} = strCellArray{2};
            loopObj = 'None';
        elseif strcmp(lineType, 'LoopProperty')
            loopPropertyNum = loopPropertyNum + 1;
            if ~strcmp(loopObj, 'Property')
                loopObj = 'Property';
                loopPropertyHead = loopPropertyNum;
            end
            crysInfo.loopProperty{loopPropertyNum, 1} = strCellArray{1};
        elseif strcmp(lineType, 'LoopValues')
            if ~strcmp(loopObj, 'Value')
                loopObj = 'Value';
                loopValueColNum = 2;
            else
                loopValueColNum = loopValueColNum + 1;
            end
            for valueIdx = 1 : numel(strCellArray)
                crysInfo.loopProperty{loopPropertyHead + valueIdx - 1, loopValueColNum} =...
                    strCellArray{valueIdx};
            end
        end
    else
        loopObj = 'None';
    end
    
    tmpTextLine = ReadTextLine;
    % in case there is one white line in the mid of the file text, if there
    % are two white lines, it is sufficient to treat it as the end of file
    if ~ischar(tmpTextLine)
        tmpTextLine = ReadTextLine;
    end
end

fclose(fileID);

% nested function:

    function [textLine] = ReadTextLine
        textLine = fgetl(fileID);
        textLine = FullString(textLine);
    end

    function [fullTextLine] = FullString(textLine)
        singleQuoteIndices = strfind(textLine, '''');
        if isempty(singleQuoteIndices)
            % no string
            fullTextLine = textLine;
        elseif length(singleQuoteIndices) == 2
            % contains full string
            fullTextLine = textLine;
        elseif length(singleQuoteIndices) == 1
            % contains incomplete string
            fullTextLine = textLine;
            nextLine = fgetl(fileID);
            singleQuoteIndices = strfind(nextLine, '''');
            while length(singleQuoteIndices) ~= 1
                fullTextLine = [fullTextLine, nextLine];
                nextLine = fgetl(fileID);
                singleQuoteIndices = strfind(nextLine, '''');
            end
            fullTextLine = [fullTextLine, nextLine];
        else
            error('Unable to parse the CIF data');
        end
    end

end



