% GnrtRandGraphene_0.m -- generate random shape graphene
clc;
close all;
clear all;
%% Generate Graphene;
a0 = 1.42;
% two distinct atoms:
vec{1} = [sqrt(3)/2, 0.5; -sqrt(3)/2, 0.5; 0, -1]';
vec{2} = [0, 1; sqrt(3)/2, -0.5; -sqrt(3) / 2, -0.5]';

GrapheneSeed = [1, 2; 0, 0; 1, 0];
IterNum = 20;
for i = 1 : IterNum
    TempGraphene = GrapheneSeed;
    for AtomIdx = 1 : size(TempGraphene, 2)
        possi = rand();
        if possi >= 0.5
            VecIdx = randperm(3);
            VecIdx = VecIdx(2);
            if TempGraphene(1, AtomIdx) == 1
                GrapheneSeed = [GrapheneSeed, [2; TempGraphene(2:3, AtomIdx) + vec{TempGraphene(1, AtomIdx)}( : ,VecIdx)]];
            else
                GrapheneSeed = [GrapheneSeed, [1; TempGraphene(2:3, AtomIdx) + vec{TempGraphene(1, AtomIdx)}( : ,VecIdx)]];
            end
        end
    end
end
GrapheneSeed = GrapheneSeed';
GrapheneSeed = uniquetol(GrapheneSeed, 1e-3, 'ByRows', true);
GrapheneSeed = GrapheneSeed';

GrapheneSeed(2 : 3, : ) = a0 * GrapheneSeed(2 : 3, : );
figure;
scatter(GrapheneSeed(2, : ), GrapheneSeed(3, : ));

Lx = max(GrapheneSeed(2, : )) - min(GrapheneSeed(2, : ));
xshift = (max(GrapheneSeed(2, : )) + min(GrapheneSeed(2, : ))) / 2;
Ly = max(GrapheneSeed(3, : )) - min(GrapheneSeed(3, : ));
yshift = (max(GrapheneSeed(3, : )) + min(GrapheneSeed(3, : ))) / 2;
L = max(Lx, Ly);
GrapheneSeed(2, : ) = GrapheneSeed(2, : ) - xshift;
GrapheneSeed(3, : ) = GrapheneSeed(3, : ) - yshift;
Lx = 1.5 * L;
Ly = 1.5 * L;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

ProjPot = ProjectedPotential(Lx, Ly, Nx, Ny, 6 * ones(1, size(GrapheneSeed, 2)), GrapheneSeed(2, : ), GrapheneSeed(3, : ));
figure;
imagesc(x, y, ProjPot);
colormap('gray'); axis square;