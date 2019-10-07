% Slicing_test_1.m -- test of RmvSlcDplAtom_0.m
clc;
close all;
clear all;

DistError = 1e-4;
%% GaAs <110>
LattConst = [3.995, 5.65];
SliceA = [31, 31, 31, 31, 33; ones(1, 5); 0, 1, 1, 0, 0.5; 0, 0, 1, 1, 0.75];
SliceB = [31, 33, 33; ones(1, 3); 0.5, 0, 1; 0.5, 0.25, 0.25];

SliceAp = RmvSlcDplAtom_0(SliceA, DistError);
SliceBp = RmvSlcDplAtom_0(SliceB, DistError);

% Comparison:
figure;
subplot(2, 2, 1);
scatter(SliceA(3, : ), SliceA(4, : ), '.r');
axis([0 1 0 1]);
subplot(2, 2, 2);
scatter(SliceAp(3, : ), SliceAp(4, : ), '.b');
axis([0 1 0 1]);
subplot(2, 2, 3);
scatter(SliceB(3, : ), SliceB(4, : ), '.r');
axis([0 1 0 1]);
subplot(2, 2, 4);
scatter(SliceBp(3, : ), SliceBp(4, : ), '.b');
axis([0 1 0 1]);

%% Test slice: non-realistic
SliceX = [25 * ones(1, 4), 30 * ones(1, 4), 45 * ones(1, 4), 57, 57, 59 * ones(1, 6);...
    0.3 * ones(1, 4), 0.2 * ones(1, 4), 0.5 * ones(1, 4), ones(1, 8);...
    0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0.25, 0.75, 0.5, 0, 1, 0, 1, 0.5;...
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0.75, 0.25, 0, 0.5, 0.5, 0.75, 0.75, 1];
SliceXp = RmvSlcDplAtom_0(SliceX, DistError);
SliceXp_1 = SliceXp( : , (SliceXp(1, : ) == 25) | (SliceXp(1, : ) == 57) | (SliceXp(1, : ) == 59));
SliceX_1 = SliceX( : , (SliceX(1, : ) == 25) | (SliceX(1, : ) == 57) | (SliceX(1, : ) == 59));
SliceXp_2 = SliceXp( : , (SliceXp(1, : ) == 30) | (SliceXp(1, : ) == 57) | (SliceXp(1, : ) == 59));
SliceX_2 = SliceX( : , (SliceX(1, : ) == 30) | (SliceX(1, : ) == 57) | (SliceX(1, : ) == 59));
SliceXp_3 = SliceXp( : , (SliceXp(1, : ) == 45) | (SliceXp(1, : ) == 57) | (SliceXp(1, : ) == 59));
SliceX_3 = SliceX( : , (SliceX(1, : ) == 45) | (SliceX(1, : ) == 57) | (SliceX(1, : ) == 59));

% comparison:
figure;
subplot(3, 2, 1);
scatter(SliceX_1(3, : ), SliceX_1(4, : ), '.r');
axis([0 1 0 1]);
subplot(3, 2, 2);
scatter(SliceXp_1(3, : ), SliceXp_1(4, : ), '.b');
axis([0 1 0 1]);
subplot(3, 2, 3);
scatter(SliceX_2(3, : ), SliceX_2(4, : ), '.r');
axis([0 1 0 1]);
subplot(3, 2, 4);
scatter(SliceXp_2(3, : ), SliceXp_2(4, : ), '.b');
axis([0 1 0 1]);
subplot(3, 2, 5);
scatter(SliceX_3(3, : ), SliceX_3(4, : ), '.r');
axis([0 1 0 1]);
subplot(3, 2, 6);
scatter(SliceXp_3(3, : ), SliceXp_3(4, : ), '.b');
axis([0 1 0 1]);