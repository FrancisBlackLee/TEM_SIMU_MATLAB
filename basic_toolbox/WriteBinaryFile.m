function WriteBinaryFile(filename, data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename, 'w');
fwrite(fileID, data, 'double');
fclose(fileID);

end

