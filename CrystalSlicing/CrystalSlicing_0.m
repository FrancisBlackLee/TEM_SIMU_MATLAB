function [Lp, SliceInfo] = CrystalSlicing_0(L, YN)
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

