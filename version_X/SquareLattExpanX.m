function [sliceOut] = SquareLattExpanX(sliceIn, lattConst, expanNum, distError, massNum, debyeTemp)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%CellNum [Nx, Ny, Nz] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   sliceIn -- fractional atomic coordinates (scaled to 0~1) and atomic 
%       types [T1, ..., TN; P1, ..., PN; x1,..., xN; y1,..., yN; z1,..., zN],
%       where T denotes atomic type, represented by the atomic numbers, P
%       denotes the elemental proportion and like before, x, y, z denote
%       the atomic coordinates;
%   lattConst -- planar lattice constants [a, b];
%   expanNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny, Nz];
%   NOTE: X denotes an experimental version!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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

if nargin == 4
    if length(expanNum) == 2
        % Transpose the matrix to improve speed
        sliceOut = sliceIn';
        coordShift = expanNum / 2 + 1; % center the lattice
        % Expand along x
        sliceBase = sliceOut;
        for i = 1: expanNum(1) + 1
            sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3) + i, sliceBase(:, 4), sliceBase(:, 5)]];
        end
        % Expand along y
        sliceBase = sliceOut;
        for i = 1: expanNum(2) + 1
            sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3), sliceBase(:, 4) + i, sliceBase(:, 5)]];
        end
        sliceOut(:, 3) = sliceOut(:, 3) - coordShift(1);
        sliceOut(:, 4) = sliceOut(:, 4) - coordShift(2);
        [row, column] = find((sliceOut(:, 3) <= expanNum(1) / 2 + distError) & (sliceOut(:, 3) >= -expanNum(1) / 2 - distError) ...
                           & (sliceOut(:, 4) <= expanNum(2) / 2 + distError) & (sliceOut(:, 4) >= -expanNum(2) / 2 - distError));
        sliceOut = sliceOut(row, :);
        sliceOut = uniquetol(sliceOut, distError, 'ByRows', true);
        sliceOut(:, 3) = lattConst(1) * sliceOut(:, 3);
        sliceOut(:, 4) = lattConst(2) * sliceOut(:, 4);
        sliceOut = sliceOut';
    else
        % Transpose the matrix to improve speed
        sliceOut = sliceIn';
        coordShift = expanNum / 2 + 1; % center the lattice
        % Expand along x
        sliceBase = sliceOut;
        for i = 1: expanNum(1) + 1
            sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3) + i, sliceBase(:, 4), sliceBase(:, 5)]];
        end
        % Expand along y
        sliceBase = sliceOut;
        for i = 1: expanNum(2) + 1
            sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3), sliceBase(:, 4) + i, sliceBase(:, 5)]];
        end
        % Expand along z
        sliceBase = sliceOut;
        for i = 1: expanNum(3) + 1
            sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3), sliceBase(:, 4), sliceBase(:, 5) + i]];
        end
        sliceOut(:, 3) = sliceOut(:, 3) - coordShift(1);
        sliceOut(:, 4) = sliceOut(:, 4) - coordShift(2);
        sliceOut(:, 5) = sliceOut(:, 5) - coordShift(3);
        [row, column] = find((sliceOut(:, 3) <= expanNum(1) / 2 + distError) & (sliceOut(:, 3) >= -expanNum(1) / 2 - distError) ...
                           & (sliceOut(:, 4) <= expanNum(2) / 2 + distError) & (sliceOut(:, 4) >= -expanNum(2) / 2 - distError) ...
                           & (sliceOut(:, 5) <= expanNum(3) / 2 + distError) & (sliceOut(:, 5) >= -expanNum(3) / 2 - distError));
        sliceOut = sliceOut(row, :);
        sliceOut = uniquetol(sliceOut, distError, 'ByRows', true);
        sliceOut(:, 3) = lattConst(1) * sliceOut(:, 3);
        sliceOut(:, 4) = lattConst(2) * sliceOut(:, 4);
        sliceOut(:, 5) = lattConst(3) * sliceOut(:, 5);
        sliceOut = sliceOut';
    end
else
    % Combine the SliceIn, MassNum and DebyeTemp:
    sliceIn = [sliceIn; massNum; debyeTemp];
    % Transpose the matrix to improve speed
    sliceOut = sliceIn';
    coordShift = expanNum / 2 + 1; % center the lattice
    distError = 1e-2;
    % Expand along x
    sliceBase = sliceOut;
    for i = 1: expanNum(1) + 1
        sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3) + i, sliceBase(:, 4), sliceBase(:, 5), sliceBase(:, 6), sliceBase(:, 7)]];
    end
    % Expand along y
    sliceBase = sliceOut;
    for i = 1: expanNum(2) + 1
        sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3), sliceBase(:, 4) + i, sliceBase(:, 5), sliceBase(:, 6), sliceBase(:, 7)]];
    end
    % Expand along z
    sliceBase = sliceOut;
    for i = 1: expanNum(3) + 1
        sliceOut = [sliceOut; [sliceBase(:, 1), sliceBase(:, 2), sliceBase(:, 3), sliceBase(:, 4), sliceBase(:, 5) + i, sliceBase(:, 6), sliceBase(:, 7)]];
    end
    sliceOut(:, 3) = sliceOut(:, 3) - coordShift(1);
    sliceOut(:, 4) = sliceOut(:, 4) - coordShift(2);
    sliceOut(:, 5) = sliceOut(:, 5) - coordShift(3);
    [row, column] = find((sliceOut(:, 3) <= expanNum(1) / 2 + distError) & (sliceOut(:, 3) >= -expanNum(1) / 2 - distError) ...
                       & (sliceOut(:, 4) <= expanNum(2) / 2 + distError) & (sliceOut(:, 4) >= -expanNum(2) / 2 - distError) ...
                       & (sliceOut(:, 5) <= expanNum(3) / 2 + distError) & (sliceOut(:, 5) >= -expanNum(3) / 2 - distError));
    sliceOut = sliceOut(row, :);
    sliceOut = uniquetol(sliceOut, distError, 'ByRows', true);
    sliceOut(:, 3) = lattConst(1) * sliceOut(:, 3);
    sliceOut(:, 4) = lattConst(2) * sliceOut(:, 4);
    sliceOut(:, 5) = lattConst(3) * sliceOut(:, 5);
    sliceOut = sliceOut';
end


end
