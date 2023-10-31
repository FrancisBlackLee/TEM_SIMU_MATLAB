% test_InitNonorthoMesh2D.m
clc;
clear;
close all;
%% main:
a = 2;
a1 = a * [1/2, sqrt(3)/2, 0.0];
a2 = a * [1/2, -sqrt(3)/2, 0.0];
na1 = 3;
na2 = 3;
n1 = 256;
n2 = 256;
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2);
rMesh = sqrt(xMesh.^2 + yMesh.^2);
frMesh = sqrt(fxMesh.^2 + fyMesh.^2);

figure;
mesh(xMesh, yMesh, rMesh);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
zlabel('r (\AA)', 'Interpreter', 'latex')

figure;
mesh(fxMesh, fyMesh, frMesh);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
zlabel('$f_r$ (1 / \AA)', 'Interpreter', 'latex')