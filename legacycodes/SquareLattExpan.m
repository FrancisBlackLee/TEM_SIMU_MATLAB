function [sliceOut] = SquareLattExpan(sliceIn, lattConst, expanNum)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%CellNum [Nx, Ny] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   SliceIn -- proportional atomic coordinates (scaled to 0~1) [x1,..., xN;
%               y1,..., yN];
%   LattConst -- planar lattice constants [a, b];
%   SliceDist -- distances between each two slices [D1,..., DN];
%   CellNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transpose the matrix to improve speed
sliceOut = sliceIn';
coordShift = expanNum / 2 + 1; % center the lattice
distError = 1e-2;
% Expand along x
sliceBase = sliceOut;
for i = 1: expanNum(1) + 1
    sliceOut = [sliceOut; [sliceBase(:, 1) + i, sliceBase(:, 2)]];
end
sliceBase = sliceOut;
for i = 1: expanNum(2) + 1
    sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2) + i]];
end
sliceOut(:, 1) = sliceOut(:, 1) - coordShift(1);
sliceOut(:, 2) = sliceOut(:, 2) - coordShift(2);
[row, ~] = find((sliceOut(:, 1) <= expanNum(1) / 2 + distError) & (sliceOut(:, 1) >= -expanNum(1) / 2 - distError) ...
                   & (sliceOut(:, 2) <= expanNum(2) / 2 + distError) & (sliceOut(:, 2) >= -expanNum(2) / 2 - distError));
sliceOut = sliceOut(row, :);
sliceOut = uniquetol(sliceOut, distError, 'ByRows', true);
sliceOut(:, 1) = lattConst(1) * sliceOut(:, 1);
sliceOut(:, 2) = lattConst(2) * sliceOut(:, 2);
sliceOut = sliceOut';

end

