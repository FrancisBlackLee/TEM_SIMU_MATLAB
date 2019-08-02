clc;
clear all;
close all;
L = load('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_010_Coordinates.txt');
L = [L(:,6),L(:,7),L(:,8),L(:,1)];
L = L';
theta_1 = -6.5; % Used to rotate the lattice
theta_1_rad = theta_1*pi/180;
theta_2 = 12.7;
theta_2_rad = theta_2*pi/180;
TransMat = eye(4);
TransMat(1:2,1:2) = [cos(theta_1_rad), -sin(theta_1_rad); sin(theta_1_rad), cos(theta_1_rad)];
L = TransMat * L;
[Lp, SliceInfo] = CrystalSlicing_0(L,0);

CrysConst = [5.43320 5.50340 23.25337];
RectPeriod = [CrysConst(1) CrysConst(3)*cos(theta_2_rad) CrysConst(2)];
CrysCenter = [(max(Lp(1,:))+min(Lp(1,:)))/2 (max(Lp(2,:))+min(Lp(2,:)))/2 (max(Lp(3,:))+min(Lp(3,:)))/2];

% Find the atoms within the rectangle denoted by RectPeriod
Rect_Index = find((Lp(1,:)>=CrysCenter(1)-RectPeriod(1)/2)&(Lp(1,:)<=CrysCenter(1)+RectPeriod(1)/2)&(Lp(2,:)>=CrysCenter(2)-RectPeriod(2)/2)&(Lp(2,:)<=CrysCenter(2)+RectPeriod(2)/2)&(Lp(3,:)>=CrysCenter(3)-RectPeriod(3)/2)&(Lp(3,:)<=CrysCenter(3)+RectPeriod(3)/2));
Unit_Rect_Latt = Lp(:,Rect_Index);
% Show the unit rectangular lattice
% Step 1: rearrange the SliceInfo array
Slice_Z = Unit_Rect_Latt(3,1);
n = 1;
Re_SliceInfo = 1;
for i = 2:size(Unit_Rect_Latt,2)
    if abs(Unit_Rect_Latt(3,i)-Slice_Z) >= 1
        Re_SliceInfo = [Re_SliceInfo 1];
        n = n+1;
        Slice_Z = Unit_Rect_Latt(3,i);
    else
        Re_SliceInfo(n) = Re_SliceInfo(n) + 1;
    end
end
% Scatter the projection of the lattice
figure;
scatter(Unit_Rect_Latt(1,:),Unit_Rect_Latt(2,:),'o','b');
axis([CrysCenter(1)-RectPeriod(1)/2 CrysCenter(1)+RectPeriod(1)/2 CrysCenter(2)-RectPeriod(2)/2 CrysCenter(2)+RectPeriod(2)/2]);
axis equal;
title('Projection of the rectangular lattice');
% Scatter the atoms on each slice
n = 1;
for i = 1:length(Re_SliceInfo)
    figure;
    hold on;
    for j = n:n+Re_SliceInfo(i)-1
        scatter(Unit_Rect_Latt(1,j),Unit_Rect_Latt(2,j),'o','b');
        text(Unit_Rect_Latt(1,j)+0.2,Unit_Rect_Latt(2,j),num2str(Unit_Rect_Latt(4,j)));
    end
    axis([CrysCenter(1)-RectPeriod(1)/2 CrysCenter(1)+RectPeriod(1)/2 CrysCenter(2)-RectPeriod(2)/2 CrysCenter(2)+RectPeriod(2)/2]);
    axis equal;
    title(['z= ' num2str(Unit_Rect_Latt(3,n))]);
    hold off;
    n = n + Re_SliceInfo(i);
end