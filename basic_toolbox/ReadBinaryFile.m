function [data] = ReadBinaryFile(filename, size)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename);
data = fread(fileID, size, 'double');
fclose(fileID);

end

