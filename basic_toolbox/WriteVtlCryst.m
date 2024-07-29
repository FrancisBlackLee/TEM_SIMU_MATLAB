function WriteVtlCryst(filename, bases, typeCoords, wobbles)
%UNTITLED2 Summary of this function goes here
%   comment
%   (ax, ay, az)
%   (bx, by, bz)
%   (cx, cy, cz)
%   type  fracX  fracY  fracZ  occu  anisoUa  anisoUb  anisoUc
%   -1

fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file!');
else
    fprintf(fid, 'Created by TEM_SIMU_MATLAB\n');
    % a
    fprintf(fid, '(%.8f  %.8f  %.8f)\n', bases.a);
    % b
    fprintf(fid, '(%.8f  %.8f  %.8f)\n', bases.b);
    % c
    fprintf(fid, '(%.8f  %.8f  %.8f)\n', bases.c);

    atomNum = size(typeCoords, 2);
    for atomIdx = 1 : atomNum
        fprintf(fid, '%d\t%.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f\n',...
            typeCoords(1, atomIdx), typeCoords(3, atomIdx), typeCoords(4, atomIdx),...
            typeCoords(5, atomIdx), typeCoords(2, atomIdx), wobbles(1, atomIdx),...
            wobbles(2, atomIdx), wobbles(3, atomIdx));
    end
    fprintf(fid, '-1\n');

    fclose(fid);
end

end