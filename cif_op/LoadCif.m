function [crysInfo] = LoadCif(filename)
% LoadCif() loads a cif file and reads the basic information from the cif
% file for other functions of TEM_SIMU_MATLAB.

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

disp('CIF content');
textLine = fgetl(fileID);
while ischar(textLine)
    [strCellArray, lineType, canDelete] = CifLineParser(textLine, loopObj);
    % disp(strCellArray);
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
    
    textLine = fgetl(fileID);
    % in case there is one white line in the mid of the file text, if there
    % are two white lines, it is sufficient to treat it as the end of file
    if ~ischar(textLine)
        textLine = fgetl(fileID);
    end
end

fclose(fileID);

end

