function [unitCell, primVecs, origin, pVecs, p] = ReadQePpOutIflag3OutputFormat3(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename);
if fid == -1
    error('Failed to open file!');
else
    % label: CRYSTAL
    line = fgetl(fid);
    % label: PRIMVEC
    line = fgetl(fid);
    % values: PRIMVEC
    primVecs = fscanf(fid, '%f', [3, 3]);
    % label: PRIMCOORD
    line = fgetl(fid);
    line = fgetl(fid);
    % nat ntyp
    nat = fscanf(fid, '%d', 1);
    unitCell = ones(5, nat);
    ntyp = fscanf(fid, '%d', 1);
    line = fgetl(fid);
    for iat = 1 : nat
        line = fgetl(fid);
        names = split(line);
        if isempty(names{1})
            names(1) = [];
        end
        unitCell(1, iat) = AtomTypeStrToIdx(names{1});
        unitCell(3, iat) = str2double(names{2});
        unitCell(4, iat) = str2double(names{3});
        unitCell(5, iat) = str2double(names{4});
    end

    % label: BEGIN_BLOCK_DATAGRID_3D
    line = fgetl(fid);
    % label: 3D_PWSCF
    line = fgetl(fid);
    % label: DATAGRID_3D_UNKNOWN
    line = fgetl(fid);
    % values: nx, ny, nz
    nx = fscanf(fid, '%d', 1);
    ny = fscanf(fid, '%d', 1);
    nz = fscanf(fid, '%d', 1);
    % values: origin
    origin = fscanf(fid, '%f', [1, 3]);
    % values: pVecs
    pVecs = fscanf(fid, '%f', [3, 3]);
    % values: p
    p = fscanf(fid, '%f', nx*ny*nz);
    p = reshape(p, [nx, ny, nz]);
    p = permute(p, [2, 1, 3]);

    fclose(fid);
end

end