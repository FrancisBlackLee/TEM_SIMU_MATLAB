%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

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
function [NewTypeCoord] = RmvSlcDplAtom_0(FracTypeCoord, DistError)
%RmvSlcDplAtom_0.m removes the periodically duplicate atoms on each slice.
%   Version 0:
%   FracTypeCoord -- fractional atomic type-coordinate matrix:
%       format: [T; P; FracX; FracY; FracZ], note that T denotes atomic
%       types represented by their atomic numbers, P denotes the atomic
%       proportions, FracX, FracY, FracZ denote the fractional atomic
%       coordinates, FracZ is not requested.

% search along x:
% sort the matrix first so that the greater duplicate coord is removed
[Xsort, Xorder] = sort(FracTypeCoord(3, : ));
NewTypeCoord = FracTypeCoord( : , Xorder);
StatIdx = 1;
while StatIdx <= size(NewTypeCoord, 2) - 1
    SrchIdx = StatIdx + 1;
    while SrchIdx <=size(NewTypeCoord, 2)
        if (abs(abs(NewTypeCoord(3, SrchIdx) - NewTypeCoord(3, StatIdx)) - 1) < DistError)&&...
                (abs(NewTypeCoord(4, SrchIdx) - NewTypeCoord(4, StatIdx)) < DistError)&&...
                (NewTypeCoord(1, SrchIdx) == NewTypeCoord(1, StatIdx))
            NewTypeCoord( : , SrchIdx) = [];
        else
            SrchIdx = SrchIdx + 1;
        end
    end
    StatIdx = StatIdx + 1;
end
% search along y:
% sort y first
[Ysort, Yorder] = sort(NewTypeCoord(4, : ));
NewTypeCoord = NewTypeCoord( : , Yorder);
StatIdx = 1;
while StatIdx <= size(NewTypeCoord, 2) - 1
    SrchIdx = StatIdx + 1;
    while SrchIdx <=size(NewTypeCoord, 2)
        if (abs(abs(NewTypeCoord(4, SrchIdx) - NewTypeCoord(4, StatIdx)) - 1) < DistError)&&...
                (abs(NewTypeCoord(3, SrchIdx) - NewTypeCoord(3, StatIdx)) < DistError)&&...
                (NewTypeCoord(1, SrchIdx) == NewTypeCoord(1, StatIdx))
            NewTypeCoord( : , SrchIdx) = [];
        else
            SrchIdx = SrchIdx + 1;
        end
    end
    StatIdx = StatIdx + 1;
end
% search along diagonal:
% sort by the distance to the origin:
[DiagSort, DiagOrder] = sort(NewTypeCoord(3, : ).^2 + NewTypeCoord(4, : ).^2);
NewTypeCoord = NewTypeCoord( : , DiagOrder);
StatIdx = 1;
while StatIdx <= size(NewTypeCoord, 2) - 1
    SrchIdx = StatIdx + 1;
    while SrchIdx <=size(NewTypeCoord, 2)
        if (abs(abs(NewTypeCoord(3, SrchIdx) - NewTypeCoord(3, StatIdx)) - 1) < DistError)&&...
           (abs(abs(NewTypeCoord(4, SrchIdx) - NewTypeCoord(4, StatIdx)) - 1) < DistError)&&...
                (NewTypeCoord(1, SrchIdx) == NewTypeCoord(1, StatIdx))
            NewTypeCoord( : , SrchIdx) = [];
        else
            SrchIdx = SrchIdx + 1;
        end
    end
    StatIdx = StatIdx + 1;
end

end

