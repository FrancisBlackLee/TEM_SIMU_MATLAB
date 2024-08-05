function [bases, typeCoords, vbases, wobbles] = ReadVtlCryst(filename)
%UNTITLED Summary of this function goes here
%   comment
%   (ax, ay, az)
%   (bx, by, bz)
%   (cx, cy, cz)
%   type  fracX  fracY  fracZ  occu  anisoUa  anisoUb  anisoUc
%   -1

fileID = fopen(filename);
if fileID == -1
    error('Failed to open file!');
else
    commentLine = fgetl(fileID);
    % a
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    ax = str2double(vecStrs(1));
    ay = str2double(vecStrs(2));
    az = str2double(vecStrs(3));
    bases.a = [ax, ay, az];

    % b
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    bx = str2double(vecStrs(1));
    by = str2double(vecStrs(2));
    bz = str2double(vecStrs(3));
    bases.b = [bx, by, bz];

    % c
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    cx = str2double(vecStrs(1));
    cy = str2double(vecStrs(2));
    cz = str2double(vecStrs(3));
    bases.c = [cx, cy, cz];

    % va
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    vax = str2double(vecStrs(1));
    vay = str2double(vecStrs(2));
    vaz = str2double(vecStrs(3));
    vbases.a = [vax, vay, vaz];

    % vb
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    vbx = str2double(vecStrs(1));
    vby = str2double(vecStrs(2));
    vbz = str2double(vecStrs(3));
    vbases.b = [vbx, vby, vbz];

    % vc
    vecLine = fgetl(fileID);
    vecLine = extractBetween(vecLine, '(', ')');
    vecStrs = split(vecLine);
    vcx = str2double(vecStrs(1));
    vcy = str2double(vecStrs(2));
    vcz = str2double(vecStrs(3));
    vbases.c = [vcx, vcy, vcz];

    tmpType = fscanf(fileID, '%d', 1);
    atomCount = 0;
    while tmpType ~= -1
        atomCount = atomCount + 1;
        typeCoords(1, atomCount) = tmpType;
        typeCoords(3, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(4, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(5, atomCount) = fscanf(fileID, '%f', 1);
        typeCoords(2, atomCount) = fscanf(fileID, '%f', 1);
        wobbles(1, atomCount) = fscanf(fileID, '%f', 1);
        wobbles(2, atomCount) = fscanf(fileID, '%f', 1);
        wobbles(3, atomCount) = fscanf(fileID, '%f', 1);
        
        tmpType = fscanf(fileID, '%d', 1);
    end

    fclose(fileID);
end

end