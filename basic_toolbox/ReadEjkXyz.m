function [lattConsts, typeCoords, wobbles] = ReadEjkXyz(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename);
if fileID == -1
    error('Failed to open file!');
else
    commentLine = fgetl(fileID);
    lattConsts(1) = fscanf(fileID, '%f', 1);
    lattConsts(2) = fscanf(fileID, '%f', 1);
    lattConsts(3) = fscanf(fileID, '%f', 1);
    tmpType = fscanf(fileID, '%d', 1);
    atomCount = 0;
    while tmpType ~= -1
        atomCount = atomCount + 1;
        typeCoords(1, atomCount) = tmpType;
        typeCoords(3, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(4, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(5, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(2, atomCount) = fscanf(fileID, '%f', 1);
        wobbles(atomCount) = fscanf(fileID, '%f', 1);
        
        tmpType = fscanf(fileID, '%d', 1);
    end
    fclose(fileID);
end

end