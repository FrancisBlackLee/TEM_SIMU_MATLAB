function [sortCrysMat, sliceInfo] = CrystalSlicing_0(crysMat, YN)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   L --input lattice data, the fourth row denotes the atom types, and the
%   first to the third row denote the atomic coordinates;
%   YN --whether to show each slice: 1 --yes, 0 --no.

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

[coordZ, zOrder] = sort(crysMat(3,:));
sortCrysMat = crysMat(:,zOrder);
sliceInfo = 1;
n = 1;
sliceZ = coordZ(1);
for i = 2:length(coordZ)
    if abs(coordZ(i)-sliceZ) >= 1e-2
        sliceInfo = [sliceInfo 1];
        n = n+1;
        sliceZ = coordZ(i);
    else
        sliceInfo(n) = sliceInfo(n) + 1;
    end
    sortCrysMat(3,i) = sliceZ;
end
% Show the slices
if YN == 1
    n = 1;
    for i = 1:length(sliceInfo)
        figure;
        hold on;
        for j = n:n+sliceInfo(i)-1
            if sortCrysMat(4,j)~=0
                scatter(sortCrysMat(1,j),sortCrysMat(2,j),'o','b');
%                 text(Lp(1,j)+0.2,Lp(2,j),num2str(Lp(4,j)));
            end
        end
        axis([min(sortCrysMat(1,:)) max(sortCrysMat(1,:)) min(sortCrysMat(2,:)) max(sortCrysMat(2,:))]);
        axis equal;
        title(['z= ' num2str(coordZ(n))]);
        n=n+sliceInfo(i);
    end
end

end

