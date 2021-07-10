function [ Lrp, SliceInfo ] = CrystalSlicing( L, t, q )
% Function CrystalSlicing slices a crystal along a certain direction
% represented by a vector t and return the results in a matrix called
% SliceInfo.
% There are also another output containing the atoms' positions after
% coordinate transformation.
% In addition, by inputting q=1 or 0 to choose whether to plot each slice
% or not.

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

% Construct a orthonormal basis containing t
A = [t(2),-t(1),0];
B = [t(1)*t(3),t(2)*t(3),-((t(1))^2+(t(2))^2)];
A = A/sqrt(sum(A.^2));
B = B/sqrt(sum(B.^2));
C = t/sqrt(sum(t.^2));
% Contruct the operator to transform coordinates
R = [A 0;B 0;C 0;0 0 0 1];
% Perform coordinate
Lr = R*L;
% Sort the matrix by its "height" in the new frame;
[Z, I] = sort(Lr(3,:));
Lrp = Lr(:,I);
SliceInfo = 1;
n = 1;
for i = 1:length(Z)-1
    if abs(Z(i)-Z(i+1)) >= 1e-2
        SliceInfo = [SliceInfo 1];
        n = n+1;
    else
        SliceInfo(n) = SliceInfo(n)+1;
    end
end
% Draw each slice
if q == 1
    n = 1;
    for i = 1:length(SliceInfo)
        figure;
        for j = n:n+SliceInfo(i)-1
            if Lrp(4,j)~=0
                hold on;
                scatter(Lrp(1,j),Lrp(2,j),'o','b');
                text(Lrp(1,j)+0.02,Lrp(2,j),num2str(Lrp(4,j)));
            end
        end
        axis([min(Lrp(1,:)) max(Lrp(1,:)) min(Lrp(2,:)) max(Lrp(2,:))]);
        axis equal;
        title(['z= ' num2str(Z(n))]);
        n=n+SliceInfo(i);
    end
end

end

