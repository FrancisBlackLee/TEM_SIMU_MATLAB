function [Lp, SliceInfo, SliceDist] = CrystalSlicing_2(L, DistError, YN)
%CrystalSlicing.m slices a given crystal described by the atomic numbers
%and atomic coordinates.
%   L -- Crystal matrix, where the first row denotes the atomic types, the
%       second row denotes the fractional concentrations and the third to
%       the fifth rows denote the atomic coordinates, whether fractional or
%       orthogonal;
%   DistError -- the largest error distance to judge whether atoms of
%       different heights be rearranged to one slice;
%   YN -- whether to show each slice: 1 --yes, 0 --no.

[Z, Order] = sort(L(5,:));
Lp = L(:,Order);
SliceInfo = 1;
n = 1;
Slice_Z = Z(1);
for i = 2:length(Z)
    if abs(Z(i)-Slice_Z) >= DistError
        SliceInfo = [SliceInfo 1];
        n = n + 1;
        Slice_Z = Z(i);
    else
        SliceInfo(n) = SliceInfo(n) + 1;
    end
    Lp(5,i) = Slice_Z;
end
SliceDist = zeros(size(SliceInfo));
n = SliceInfo(1);
for i = 1 : length(SliceInfo) - 1
    SliceDist(i) = Lp(5, n + 1) - Lp(5, n);
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
            if Lp(1,j)~=0
                scatter(Lp(3,j),Lp(4,j),'o','b');
%                 text(Lp(1,j)+0.2,Lp(2,j),num2str(Lp(4,j)));
            end
        end
        axis([min(Lp(3,:)) max(Lp(3,:)) min(Lp(4,:)) max(Lp(4,:))]);
        axis equal;
        title(['z= ' num2str(Z(n))]);
        n=n+SliceInfo(i);
    end
end

end

