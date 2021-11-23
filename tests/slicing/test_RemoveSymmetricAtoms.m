% test_RemoveSymmetricAtoms.m
clc;
clear;
close all;
%% test 1:
fracCoords_1 = [6, 1.0000, 0.0000, 0.0000, 0.0000;...
              6, 1.0000, 1.0000, 1.0000, 1.0000;...
              6, 1.0000, 0.5000, 0.5000, 0.5000;...
              6, 1.0000, 0.5000, 0.5000, 0.0000;...
              6, 1.0000, 0.5000, 0.5000, 1.0000]';
figure;
subplot(1, 2, 1)
scatter3(fracCoords_1(3, :), fracCoords_1(4, :), fracCoords_1(5, :), 'filled');
title('test 1: old');
xlabel('x');
xlim([0, 1]);
ylabel('y');
ylim([0, 1]);
zlabel('z');
zlim([0, 1]);

tolerance = 1.0e-8;
newCoords_1 = RemoveSymmetricAtoms(fracCoords_1, tolerance);
subplot(1, 2, 2)
scatter3(newCoords_1(3, :), newCoords_1(4, :), newCoords_1(5, :), 'filled');
title('test 1: new');
xlabel('x');
xlim([0, 1]);
ylabel('y');
ylim([0, 1]);
zlabel('z');
zlim([0, 1]);

%% test 2:
fracCoords_2 = [6, 0.6000, 0.0000, 0.0000, 0.0000;...
                6, 0.6000, 1.0000, 1.0000, 1.0000;...
                6, 0.6000, 0.5000, 0.5000, 0.5000;...
                6, 0.6000, 0.5000, 0.5000, 0.0000;...
                6, 0.6000, 0.5000, 0.5000, 1.0000;...
                7, 0.4000, 0.0000, 0.0000, 0.0000;...
                7, 0.4000, 1.0000, 1.0000, 1.0000;...
                7, 0.4000, 0.5000, 0.5000, 0.5000;...
                7, 0.4000, 0.5000, 0.5000, 0.0000;...
                7, 0.4000, 0.5000, 0.5000, 1.0000]';
figure;
subplot(1, 2, 1)
scatter3(fracCoords_2(3, :), fracCoords_2(4, :), fracCoords_2(5, :), 'filled');
title('test 2: old');
xlabel('x');
xlim([0, 1]);
ylabel('y');
ylim([0, 1]);
zlabel('z');
zlim([0, 1]);

tolerance = 1.0e-8;
newCoords_2 = RemoveSymmetricAtoms(fracCoords_2, tolerance);
subplot(1, 2, 2)
scatter3(newCoords_2(3, :), newCoords_2(4, :), newCoords_2(5, :), 'filled');
title('test 2: new');
xlabel('x');
xlim([0, 1]);
ylabel('y');
ylim([0, 1]);
zlabel('z');
zlim([0, 1]);