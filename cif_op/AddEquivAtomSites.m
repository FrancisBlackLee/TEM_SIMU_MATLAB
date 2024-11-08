function [newAtomSiteMat] = AddEquivAtomSites(atomSiteMat, nCell, tolerance)
%AddEquivAtomSites() adds equivalent atom sites on the cell boundaries.
% Input:
%   atomSiteMat -- atomic site matrix to which equivalent atoms are added,
%       format: [type; occupancy; fractX; fractY; fractZ];
% Output:
%   newAtomSiteMat -- atomic site matrix to which equivalent atoms are
%       already added, format: [type; occupancy; fractX; fractY; fractZ];
%
% NOTE:
%   The input cell must be a unit cell.

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

if nargin == 1
    nCell = [1, 1, 1];
    tolerance = 1e-8;
elseif nargin == 2
    tolerance = 1e-8;
end

% initAtomNum = size(atomSiteMat, 2);
% 
% newAtomSiteMat = repmat(atomSiteMat, 1, 27);
% for xGlide = -1 : 1
%     start_x = (xGlide + 1) * initAtomNum * 9;
%     for yGlide = -1 : 1
%         start_y = (yGlide + 1) * initAtomNum * 3;
%         for zGlide = -1 : 1
%             start_z = (zGlide + 1) * initAtomNum;
%             startIdx = start_x + start_y + start_z + 1;
%             range = startIdx : startIdx + initAtomNum - 1;
%             newAtomSiteMat(3, range) = newAtomSiteMat(3, range) + xGlide;
%             newAtomSiteMat(4, range) = newAtomSiteMat(4, range) + yGlide;
%             newAtomSiteMat(5, range) = newAtomSiteMat(5, range) + zGlide;
%         end
%     end
% end

% x boundary
xBdAtoms = atomSiteMat(:, mod(atomSiteMat(3, :), nCell(1)) < tolerance);
xBdAtoms(3, :) = nCell(1) - xBdAtoms(3, :);

newAtomSiteMat = [atomSiteMat, xBdAtoms];

% y boundary
yBdAtoms = newAtomSiteMat(:, mod(newAtomSiteMat(4, :), nCell(2)) < tolerance);
yBdAtoms(4, :) = nCell(2) - yBdAtoms(4, :);

newAtomSiteMat = [newAtomSiteMat, yBdAtoms];

% z boundary
zBdAtoms = newAtomSiteMat(:, mod(newAtomSiteMat(5, :), nCell(3)) < tolerance);
zBdAtoms(5, :) = nCell(3) - zBdAtoms(5, :);

newAtomSiteMat = [newAtomSiteMat, zBdAtoms];

% index = find(newAtomSiteMat(3, :) > 1 + tolerance);
% newAtomSiteMat(3, index) = mod(newAtomSiteMat(3, index), 1);
% index = find(newAtomSiteMat(3, :) < 0 - tolerance);
% newAtomSiteMat(3, index) = newAtomSiteMat(3, index) - floor(newAtomSiteMat(3, index));
% 
% index = find(newAtomSiteMat(4, :) > 1 + tolerance);
% newAtomSiteMat(4, index) = mod(newAtomSiteMat(4, index), 1);
% index = find(newAtomSiteMat(4, :) < 0 - tolerance);
% newAtomSiteMat(4, index) = newAtomSiteMat(4, index) - floor(newAtomSiteMat(4, index));
% 
% index = find(newAtomSiteMat(5, :) > 1 + tolerance);
% newAtomSiteMat(5, index) = mod(newAtomSiteMat(5, index), 1);
% index = find(newAtomSiteMat(5, :) < 0 - tolerance);
% newAtomSiteMat(5, index) = newAtomSiteMat(5, index) - floor(newAtomSiteMat(5, index));

newAtomSiteMat = newAtomSiteMat';
newAtomSiteMat = uniquetol(newAtomSiteMat, tolerance, 'ByRows', true);
newAtomSiteMat = newAtomSiteMat';

end



