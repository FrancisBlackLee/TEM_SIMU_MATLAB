clc;
clear all;
close all;
CrysConst = [5.43320 5.50340 23.25337];
%% [010]
L_010 = load('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_010_Coordinates.txt');
L_010 = [L_010(:,6),L_010(:,7),L_010(:,8),L_010(:,1)];
L_010 = L_010';
theta_010_1 = -6.5; % Used to rotate the lattice
theta_010_1_rad = theta_010_1*pi/180;
theta_010_2 = 12.7;
theta_010_2_rad = theta_010_2*pi/180;
TransMat_010 = eye(4);
TransMat_010(1:2,1:2) = [cos(theta_010_1_rad), -sin(theta_010_1_rad); sin(theta_010_1_rad), cos(theta_010_1_rad)];
L_010 = TransMat_010 * L_010;

CrysCenter_010 = [(max(L_010(1,:))+min(L_010(1,:)))/2, (max(L_010(2,:))+min(L_010(2,:)))/2, (max(L_010(3,:))+min(L_010(3,:)))/2];
RectPeriod_010 = [CrysConst(1), CrysConst(3)*cos(theta_010_2_rad), CrysConst(2)];
% Find the atoms within the rectangle denoted by RectPeriod
Rect_Index_010 = find((L_010(1,:)>=CrysCenter_010(1)-RectPeriod_010(1)/2)...
                     &(L_010(1,:)<=CrysCenter_010(1)+RectPeriod_010(1)/2)...
                     &(L_010(2,:)>=CrysCenter_010(2)-RectPeriod_010(2)/2)...
                     &(L_010(2,:)<=CrysCenter_010(2)+RectPeriod_010(2)/2)...
                     &(L_010(3,:)>=CrysCenter_010(3)-RectPeriod_010(3)/2)...
                     &(L_010(3,:)<=CrysCenter_010(3)+RectPeriod_010(3)/2));
Unit_Rect_Latt_010 = L_010(:,Rect_Index_010);
% Slice the rectangular unit cell
[Unit_010, SliceInfo_010] = CrystalSlicing_0(Unit_Rect_Latt_010, 0);
Unit_Prop_010 = Unit_010;
Unit_Prop_010(1,:) = Unit_010(1,:)/RectPeriod_010(1) - CrysCenter_010(1)/RectPeriod_010(1)*ones(size(Unit_010(1,:))) + 0.5*ones(size(Unit_010(1,:)));
Unit_Prop_010(2,:) = Unit_010(2,:)/RectPeriod_010(2) - CrysCenter_010(2)/RectPeriod_010(2)*ones(size(Unit_010(1,:))) + 0.5*ones(size(Unit_010(1,:)));
Unit_010 = Unit_010';
Unit_Prop_010 = Unit_Prop_010';
SliceInfo_010 = SliceInfo_010';
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_010\LPCMO_Unit_010.txt','Unit_010','-ascii','-tabs');
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_010\LPCMO_Unit_Prop_010.txt','Unit_Prop_010','-ascii','-tabs');
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_010\LPCMO_SliceInfo_010.txt','SliceInfo_010','-ascii','-tabs');
%% [110]
L_110 = load('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_110_Coordinates.txt');
L_110 = [L_110(:,6),L_110(:,7),L_110(:,8),L_110(:,1)];
L_110 = L_110';
theta_110_1 = -9.2;
theta_110_1_rad = theta_110_1*pi/180;
TransMat_110 = eye(4);
TransMat_110(1:2,1:2) = [cos(theta_110_1_rad), -sin(theta_110_1_rad); sin(theta_110_1_rad), cos(theta_110_1_rad)];
L_110 = TransMat_110 * L_110;

CrysCenter_110 = [(max(L_110(1,:))+min(L_110(1,:)))/2, (max(L_110(2,:))+min(L_110(2,:)))/2, (max(L_110(3,:))+min(L_110(3,:)))/2];
RectPeriod_110 = [CrysConst(1)*CrysConst(2)/sqrt(CrysConst(1)^2+CrysConst(2)^2), CrysConst(3)*cos(theta_110_1_rad), sqrt(CrysConst(1)^2+CrysConst(2)^2)];
Rect_Index_110 = find((L_110(1,:)>=CrysCenter_110(1)-RectPeriod_110(1)/2)...
                     &(L_110(1,:)<=CrysCenter_110(1)+RectPeriod_110(1)/2)...
                     &(L_110(2,:)>=CrysCenter_110(2)-RectPeriod_110(2)/2)...
                     &(L_110(2,:)<=CrysCenter_110(2)+RectPeriod_110(2)/2)...
                     &(L_110(3,:)>=CrysCenter_110(3)-RectPeriod_110(3)/2)...
                     &(L_110(3,:)<=CrysCenter_110(3)+RectPeriod_110(3)/2));
Unit_Rect_Latt_110 = L_110(:,Rect_Index_110);
% Slicing the rectangular unit cell
[Unit_110, SliceInfo_110] = CrystalSlicing_0(Unit_Rect_Latt_110, 0);
Unit_Prop_110 = Unit_110;
Unit_Prop_110(1,:) = Unit_110(1,:)/RectPeriod_110(1) - CrysCenter_110(1)/RectPeriod_110(1)*ones(size(Unit_110(1,:))) + 0.5*ones(size(Unit_110(1,:)));
Unit_Prop_110(2,:) = Unit_110(2,:)/RectPeriod_110(2) - CrysCenter_110(2)/RectPeriod_110(2)*ones(size(Unit_110(1,:))) + 0.5*ones(size(Unit_110(1,:)));
Unit_110 = Unit_110';
Unit_Prop_110 = Unit_Prop_110';
SliceInfo_110 = SliceInfo_110';
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_110\LPCMO_Unit_110.txt','Unit_110','-ascii','-tabs');
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_110\LPCMO_Unit_Prop_110.txt','Unit_Prop_110','-ascii','-tabs');
% save('D:\Francis. B. Lee\cooperation\Group Cooperation\P_CLZ\LPCMO_110\LPCMO_SliceInfo_110.txt','SliceInfo_110','-ascii','-tabs');