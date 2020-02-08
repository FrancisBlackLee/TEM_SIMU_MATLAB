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
function [Lp, SliceInfo, SliceDist] = CrystalSlicing_1(L, YN)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   L --input lattice data, the fourth row denotes the atom types, and the
%   first to the third row denote the atomic coordinates;
%   YN --whether to show each slice: 1 --yes, 0 --no.

[Z, Order] = sort(L(3,:));
Lp = L(:,Order);
SliceInfo = 1;
n = 1;
Slice_Z = Z(1);
for i = 2:length(Z)
    if abs(Z(i)-Slice_Z) >= 1e-2
        SliceInfo = [SliceInfo 1];
        n = n+1;
        Slice_Z = Z(i);
    else
        SliceInfo(n) = SliceInfo(n) + 1;
    end
    Lp(3,i) = Slice_Z;
end
SliceDist = zeros(size(SliceInfo));
n = SliceInfo(1);
for i = 1 : length(SliceInfo) - 1
    SliceDist(i) = Lp(3, n + 1) - Lp(3, n);
    n = n + SliceInfo(i + 1);
end
SliceDist(i + 1) = SliceDist(1);
% Show the slices
if YN == 1
    n = 1;
    for i = 1:length(SliceInfo)
        figure;
        hold on;
        for j = n:n+SliceInfo(i)-1
            if Lp(4,j)~=0
                scatter(Lp(1,j),Lp(2,j),'o','b');
%                 text(Lp(1,j)+0.2,Lp(2,j),num2str(Lp(4,j)));
            end
        end
        axis([min(Lp(1,:)) max(Lp(1,:)) min(Lp(2,:)) max(Lp(2,:))]);
        axis equal;
        title(['z= ' num2str(Z(n))]);
        n=n+SliceInfo(i);
    end
end

end

