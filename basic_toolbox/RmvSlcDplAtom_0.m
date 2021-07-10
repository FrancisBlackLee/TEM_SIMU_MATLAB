function [newTypeCoord] = RmvSlcDplAtom_0(fracTypeCoord, tolerance)
%RmvSlcDplAtom_0.m removes the periodically duplicate atoms on each slice.
%   Version 0:
%   fracTypeCoord -- fractional atomic type-coordinate matrix:
%       format: [T; P; FracX; FracY; FracZ], note that T denotes atomic
%       types represented by their atomic numbers, P denotes the atomic
%       proportions, FracX, FracY, FracZ denote the fractional atomic
%       coordinates, FracZ is not requested.

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

% search along x:
% sort the matrix first so that the greater duplicate coord is removed
[~, xOrder] = sort(fracTypeCoord(3, : ));
newTypeCoord = fracTypeCoord( : , xOrder);
statIdx = 1;
while statIdx <= size(newTypeCoord, 2) - 1
    srchIdx = statIdx + 1;
    while srchIdx <=size(newTypeCoord, 2)
        if (abs(abs(newTypeCoord(3, srchIdx) - newTypeCoord(3, statIdx)) - 1) < tolerance)&&...
                (abs(newTypeCoord(4, srchIdx) - newTypeCoord(4, statIdx)) < tolerance)&&...
                (newTypeCoord(1, srchIdx) == newTypeCoord(1, statIdx))
            newTypeCoord( : , srchIdx) = [];
        else
            srchIdx = srchIdx + 1;
        end
    end
    statIdx = statIdx + 1;
end
% search along y:
% sort y first
[~, yOrder] = sort(newTypeCoord(4, : ));
newTypeCoord = newTypeCoord( : , yOrder);
statIdx = 1;
while statIdx <= size(newTypeCoord, 2) - 1
    srchIdx = statIdx + 1;
    while srchIdx <=size(newTypeCoord, 2)
        if (abs(abs(newTypeCoord(4, srchIdx) - newTypeCoord(4, statIdx)) - 1) < tolerance)&&...
                (abs(newTypeCoord(3, srchIdx) - newTypeCoord(3, statIdx)) < tolerance)&&...
                (newTypeCoord(1, srchIdx) == newTypeCoord(1, statIdx))
            newTypeCoord( : , srchIdx) = [];
        else
            srchIdx = srchIdx + 1;
        end
    end
    statIdx = statIdx + 1;
end
% search along diagonal:
% sort by the distance to the origin:
[~, diagOrder] = sort(newTypeCoord(3, : ).^2 + newTypeCoord(4, : ).^2);
newTypeCoord = newTypeCoord( : , diagOrder);
statIdx = 1;
while statIdx <= size(newTypeCoord, 2) - 1
    srchIdx = statIdx + 1;
    while srchIdx <=size(newTypeCoord, 2)
        if (abs(abs(newTypeCoord(3, srchIdx) - newTypeCoord(3, statIdx)) - 1) < tolerance)&&...
           (abs(abs(newTypeCoord(4, srchIdx) - newTypeCoord(4, statIdx)) - 1) < tolerance)&&...
                (newTypeCoord(1, srchIdx) == newTypeCoord(1, statIdx))
            newTypeCoord( : , srchIdx) = [];
        else
            srchIdx = srchIdx + 1;
        end
    end
    statIdx = statIdx + 1;
end

end

