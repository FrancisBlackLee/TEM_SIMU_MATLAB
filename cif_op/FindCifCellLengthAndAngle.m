function [cellLength, cellAngle] = FindCifCellLengthAndAngle(filename)
%FindCifCellLengthAndAngle() finds cell lengths and angles from the CIF
%file directly.
% Input:
%   filename -- CIF filename;
% Output:
%   cellLengths -- [a, b, c];
%   cellAngles -- [alpha, beta, gamma];

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

