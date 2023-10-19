function [empadData] = ReadEmpadData(filename, tNum, zNum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[fileID, errMsg] = fopen(filename);
if fileID == -1
    error(errMsg);
end

imageSize = [128, 130];
cropRangeY = 3 : 126;
cropRangeX = 3 : 126;

empadData = zeros(124, 124, tNum, zNum);
for tIdx = 1 : tNum
    for zIdx = 1 : zNum
        tmpRawImage = fread(fileID, imageSize, 'float', 'ieee-le');
        empadData(:, :, tIdx, zIdx) = tmpRawImage(cropRangeY, cropRangeX)';
    end
end

fclose(fileID);

end