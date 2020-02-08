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
function [slice, SliceDist] = CrystalSlicing_3(L, DistError, Zmax, YN, PlotColor)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   L -- Crystal matrix, where the first row denotes the atomic types, the
%       second row denotes the fractional concentrations and the third to
%       the fifth rows denote the atomic coordinates, whether fractional or
%       orthogonal;
%   DistError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.
%   NOTE: this version was derived from the previous experimental version,
%   designed for LPCMO slicing.

AtomNum = size(L, 2);
[Z, Order] = sort(L(5,:));
Lp = L(:,Order);
PlotColor = PlotColor(Order);
SliceInfo = 1;
n = 1;
slice{n} = Lp( : , 1);
Slice_Z = Z(1);
for i = 2:length(Z)
    if abs(Z(i)-Slice_Z) >= DistError
        SliceInfo = [SliceInfo 1];
        n = n + 1;
        slice{n} = Lp( : , i);
        Slice_Z = Z(i);
    else
        SliceInfo(n) = SliceInfo(n) + 1;
        Lp(5, i) = Slice_Z;
        slice{n} = [slice{n}, Lp( : , i)];
    end
end
SliceDist = zeros(size(SliceInfo));
n = SliceInfo(1);
for i = 1 : length(SliceInfo) - 1
    SliceDist(i) = Lp(5, n + 1) - Lp(5, n);
    n = n + SliceInfo(i + 1);
end
if Zmax - sum(SliceDist(1 : i)) >= DistError
    SliceDist(i + 1) = Zmax - sum(SliceDist(1 : i));
else
    SliceDist(i) = SliceDist(i) + Zmax - sum(SliceDist(1 : i));
    slice{i+1}(5, : ) = Z(1);
    Lp(5, AtomNum - SliceInfo(i + 1) + 1 : AtomNum) = Lp(5, 1);
    Lp = [Lp( : , AtomNum - SliceInfo(i + 1) + 1 : AtomNum), Lp( : , 1 : AtomNum - SliceInfo(i + 1))];
    PlotColor = [PlotColor(AtomNum - SliceInfo(i + 1) + 1 : AtomNum), PlotColor(1 : AtomNum - SliceInfo(i + 1))];
    slice{1} = [slice{1}, slice{i+1}];
    slice(i + 1) = [];
    SliceDist(i + 1) = [];
    SliceInfo(1) = SliceInfo(1) + SliceInfo(i + 1);
    SliceInfo(i + 1) = [];
end
% Show the slices
if YN == 1
    n = 1;
    for i = 1:length(SliceInfo)
        figure;
        hold on;
        for j = n:n+SliceInfo(i)-1
            if Lp(1,j)~=0
                % No more than 8 types of color
                Colors = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k'];
                scatter(Lp(3, j), Lp(4, j), 'o', Colors(PlotColor(j)));
            end
        end
        axis([min(Lp(3,:)) max(Lp(3,:)) min(Lp(4,:)) max(Lp(4,:))]);
        axis equal;
        title(['z= ' num2str(Z(n))]);
        n=n+SliceInfo(i);
    end
end

end

