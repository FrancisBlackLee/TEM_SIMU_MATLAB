function [SliceOut] = SquareLattExpan_1(SliceIn, LattConst, CellNum)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%CellNum [Nx, Ny, Nz] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   SliceIn -- proportional atomic coordinates (scaled to 0~1) and atomic 
%   types [x1,..., xN; y1,..., yN; z1,..., zN; T1,..., TN];
%   LattConst -- planar lattice constants [a, b];
%   SliceDist -- distances between each two slices [D1,..., DN];
%   CellNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny, Nz];

% Transpose the matrix to improve speed
SliceOut = SliceIn';
CoordShift = CellNum / 2 + 1; % center the lattice
DistError = 1e-2;
% Expand along x
SliceBase = SliceOut;
for i = 1: CellNum(1) + 1
    SliceOut = [SliceOut; [SliceBase(:, 1) + i, SliceBase(:, 2), SliceBase(:, 3), SliceBase(:, 4)]];
end
% Expand along y
SliceBase = SliceOut;
for i = 1: CellNum(2) + 1
    SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2) + i, SliceBase(:, 3), SliceBase(:, 4)]];
end
% Expand along z
SliceBase = SliceOut;
for i = 1: CellNum(3) + 1
    SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3) + i, SliceBase(:, 4)]];
end
SliceOut(:, 1) = SliceOut(:, 1) - CoordShift(1);
SliceOut(:, 2) = SliceOut(:, 2) - CoordShift(2);
SliceOut(:, 3) = SliceOut(:, 3) - CoordShift(3);
[row, column] = find((SliceOut(:, 1) <= CellNum(1) / 2 + DistError) & (SliceOut(:, 1) >= -CellNum(1) / 2 - DistError) ...
                   & (SliceOut(:, 2) <= CellNum(2) / 2 + DistError) & (SliceOut(:, 2) >= -CellNum(2) / 2 - DistError) ...
                   & (SliceOut(:, 3) <= CellNum(3) / 2 + DistError) & (SliceOut(:, 3) >= -CellNum(3) / 2 - DistError));
SliceOut = SliceOut(row, :);
SliceOut = uniquetol(SliceOut, DistError, 'ByRows', true);
SliceOut(:, 1) = LattConst(1) * SliceOut(:, 1);
SliceOut(:, 2) = LattConst(2) * SliceOut(:, 2);
SliceOut(:, 3) = LattConst(3) * SliceOut(:, 3);
SliceOut = SliceOut';

end
