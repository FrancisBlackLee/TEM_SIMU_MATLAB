function [str] = ReadStrFromCif(fileID)
%ReadStrFromCif() read a valid string from the CIF file.
% Input:
%   fileID -- file identity/handle/pointer of the file;
% Output:
%   str -- string;

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

