function WriteVtlXyz(filename, lattConsts, typeCoords, wobbles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file!');
else
    fprintf(fid, 'Created by TEM_SIMU_MATLAB\n');
    fprintf(fid, '\t%.6f\t%.6f\t%.6f\n', lattConsts(1), lattConsts(2),...
        lattConsts(3));
    atomNum = size(typeCoords, 2);
    for atomIdx = 1 : atomNum
        fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',...
            typeCoords(1, atomIdx), typeCoords(3, atomIdx), typeCoords(4, atomIdx),...
            typeCoords(5, atomIdx), typeCoords(2, atomIdx), wobbles(1, atomIdx),...
            wobbles(2, atomIdx), wobbles(3, atomIdx));
    end
    fprintf(fid, '-1\n');

    fclose(fid);
end

end