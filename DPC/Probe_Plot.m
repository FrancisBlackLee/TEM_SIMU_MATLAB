clc
close all
clear all

L = 10;
M = 512;
dx = L/(M-1);
x = linspace(-L/2,L/2,M);
y = x;

KeV = 200;
lambda = 12.3986/sqrt((2*511.0+KeV)*KeV);
k = 2*pi/lambda;
% probe1 params:
Params1.KeV = KeV;
Params1.amax = 10.37;
Params1.Cs = 1;
Params1.df = 750;
% probe2 params:
Params2.KeV = KeV;
Params2.amax = 10.37;
Params2.Cs = 1;
Params2.df = 500;
% probe1 params:
Params3.KeV = KeV;
Params3.amax = 10.37;
Params3.Cs = 1;
Params3.df = 250;

[X,Y] = meshgrid(x,y);
fx = linspace(-1/(2*dx),1/(2*dx),M);
fy = fx;
[Fx, Fy] = meshgrid(fx,fy);

% Draw probe
Probe_test = [ProbeCreate(Params1,0,0,L,L,M,M), ProbeCreate(Params2,0,0,L,L,M,M), ProbeCreate(Params3,0,0,L,L,M,M)];
Probe_test_Intensity = abs(Probe_test.^2);
figure;
mesh(Probe_test_Intensity);
axis off;