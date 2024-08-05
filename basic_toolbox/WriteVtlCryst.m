function WriteVtlCryst(filename, bases, typeCoords, vbases, wobbles)
%UNTITLED2 Summary of this function goes here
%   comment
%   a = (ax, ay, az)
%   b = (bx, by, bz)
%   c = (cx, cy, cz)
%   va = (vax, vay, vaz)
%   vb = (vbx, vby, vbz)
%   vc = (vcx, vcy, vcz)
%   type  fracX  fracY  fracZ  occu  anisoUa  anisoUb  anisoUc
%   -1

fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file!');
else
    fprintf(fid, 'Created by TEM_SIMU_MATLAB\n');
    % a
    fprintf(fid, 'a = (%.8f  %.8f  %.8f)\n', bases.a);
    % b
    fprintf(fid, 'b = (%.8f  %.8f  %.8f)\n', bases.b);
    % c
    fprintf(fid, 'c = (%.8f  %.8f  %.8f)\n', bases.c);
    % va
    fprintf(fid, 'va = (%.8f  %.8f  %.8f)\n', vbases.a);
    % vb
    fprintf(fid, 'vb = (%.8f  %.8f  %.8f)\n', vbases.b);
    % vc
    fprintf(fid, 'vc = (%.8f  %.8f  %.8f)\n', vbases.c);

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