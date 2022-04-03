% si_111_cbed_data_view.m
clc;
clear;
close all;
%% load data:
Nx = 512;
Ny = 512;
cbeds = VtemlabImageReader('test_cbed.txt', Ny, Nx);

%% view data:
cbedIdx = 10;
rawCbed = cbeds(:, :, cbedIdx);
logCoeff = 300;
logCbed = log(1 + logCoeff * rawCbed / max(rawCbed, [], 'all'));

figure;
subplot(1, 2, 1);
imagesc(rawCbed);
axis square;

subplot(1, 2, 2);
imagesc(logCbed);
axis square;