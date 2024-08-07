function [sliceOut] = SquareLattExpan_0(sliceIn, lattConst, expanNum)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%expanNum [Nx, Ny] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   sliceIn -- proportional atomic coordinates (scaled to 0~1) and atomic 
%       types [T1,..., TN; x1,..., xN; y1,..., yN; z1,..., zN];
%   lattConst -- planar lattice constants [a, b];
%   sliceDist -- distances between each two slices [D1,..., DN];
%   expanNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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
SliceBase = sliceOut;
for i = 1: expanNum(1) + 1
    sliceOut = [sliceOut; [SliceBase(:, 1), SliceBase(:, 2) + i, SliceBase(:, 3), SliceBase(:, 4)]];
end
SliceBase = sliceOut;
for i = 1: expanNum(2) + 1
    sliceOut = [sliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3) + i, SliceBase(:, 4)]];
end
sliceOut(:, 2) = sliceOut(:, 2) - coordShift(1);
sliceOut(:, 3) = sliceOut(:, 3) - coordShift(2);
[row, ~] = find((sliceOut(:, 2) <= expanNum(1) / 2 + distError) & (sliceOut(:, 2) >= -expanNum(1) / 2 - distError) ...
                   & (sliceOut(:, 3) <= expanNum(2) / 2 + distError) & (sliceOut(:, 3) >= -expanNum(2) / 2 - distError));
sliceOut = sliceOut(row, :);
sliceOut = uniquetol(sliceOut, distError, 'ByRows', true);
sliceOut(:, 2) = lattConst(1) * sliceOut(:, 2);
sliceOut(:, 3) = lattConst(2) * sliceOut(:, 3);
sliceOut(:, 4) = lattConst(3) * sliceOut(:, 4);
sliceOut = sliceOut';

end



