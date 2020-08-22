%% Slicing_goldNC_0.m 2020/04/05
% 4nm gold in water
clc;
close all;
clear all;
%% main:
latt = load('GoldClusterCoord_4nm.txt');
latt = [latt(:, 1), ones(size(latt(:, 1))), latt(:, 6), latt(:, 7), latt(:, 8)]';

% figure;
% scatter3(latt(3, : ), latt(4, : ), latt(5, : ));

% move the nano cluster to the center of vacuum box
latt(3:4, :) = latt(3:4, :) + 30.0;
latt(5, : ) = latt(5, : ) + 50.0;
% figure;
% scatter3(latt(3, : ), latt(4, : ), latt(5, : ));


water = load('E:\Group_Affairs\Zhang_Jiarui\water_20nm.txt');
water = water';
waterAN = size(water, 2);
water = [water(1, :); ones(1, waterAN); water(2:4, :)];
% figure;
% scatter3(water(3, :), water(4, :), water(5, :));

water(:, (water(3, :)>60) | (water(4, :)>60)) = [];
water(:, water(5, :)>100) = [];
water(:, sqrt((water(3,:)-30).^2 + (water(4,:)-30).^2 + (water(5,:)-50).^2) < 20) = [];
% figure;
% scatter3(water(3, :), water(4, :), water(5, :));

sample = [latt, water];

[slice, SliceDist, ExtraSlice] = CrystalSlicing_X(sample, sample, 2.0, 100.0, 1, 0);
save('Tests\SliceDist_4nmGold_in_water.txt', 'SliceDist',  '-ascii', '-double', '-tabs');

Lx = 60.0;
Ly = Lx;
Nx = 1024;
Ny = 1024;
CellNum = [1, 1];
LattConst = [Lx, Ly];
for i = 1:length(SliceDist)
    filename = strcat('Tests\ProjPot\p', num2str(i), '.txt');
    tempSlice = slice{1, i};
    tempSlice(3, : ) = tempSlice(3, : ) / Lx;
    tempSlice(4, : ) = tempSlice(4, : ) / Ly;
    ProjPot = MultiProjPot_conv_X(tempSlice, CellNum, LattConst, Lx, Ly, Nx, Ny, 1.0e-5);
    save(filename, 'ProjPot', '-ascii', '-double', '-tabs');
    disp(i);
end