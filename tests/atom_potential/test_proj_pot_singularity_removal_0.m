% test_proj_pot_singularity_removal_0.m
clc;
clear;
close all;
%% main:
dx = 0.01;
Nx = 256;
Lx = dx * Nx;
x = -Lx / 2 : dx : Lx / 2 - dx;

frac = 1;
threshold = frac * dx;

atomType = 1;
atomNum = 20;
rndX = 0.9 * (Lx * rand(1, atomNum) - Lx / 2);
atomY = 0;
figure;
hold on;
for atomIdx = 1 : atomNum
    tmpProjPot = ProjPotTestKernel(Lx, Lx, Nx, Nx, atomType, rndX(atomIdx),...
        atomY, threshold);
    plot(x, tmpProjPot(Nx / 2 + 1, :));
end
hold off;
box on;
xlabel('$x (\AA)$', 'interpreter', 'latex');
ylabel('$V_z(x) (V)$', 'interpreter', 'latex');

folder = strcat('tests\atom_potential\threshold_data\Z=', num2str(atomType));
mkDirStat = mkdir(folder);
filename = strcat('dx=', num2str(dx, '%.2f'), '_frac=', num2str(frac, '%.2f'),...
    '.jpg');

%% extra functions:
function projPot = ProjPotTestKernel(Lx, Ly, Nx, Ny, atomType,...
    atomX, atomY, threshold)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomZ -- atomic numbers of the input atom series;
%   AtomX, AtomY -- atomic coordinates corresponding to the atomic series;

atomNum = length(atomType);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);

bohrRadius = 0.529; % Bohr radius in angstrom
eleCharge = 14.4; % elemental charge in volt - angstrom

scattParam = load('EJKScattParam.txt');

projPot = zeros(size(X));
deltaRhoSq = threshold^2;
for atomIdx = 1 : atomNum
    rhoSq = (X - atomX(atomIdx)).^2 + (Y - atomY(atomIdx)).^2;
    rhoSq(rhoSq < deltaRhoSq) = deltaRhoSq;
    startIndex = 3 * (atomType(atomIdx) - 1) + 1;
    paramA = [scattParam(startIndex, 1),...
        scattParam(startIndex, 3),...
        scattParam(startIndex + 1, 1)];
    paramB = [scattParam(startIndex, 2),...
        scattParam(startIndex, 4),...
        scattParam(startIndex + 1, 2)];
    paramC = [scattParam(startIndex + 1, 3),...
        scattParam(startIndex + 2, 1),...
        scattParam(startIndex + 2, 3)];
    paramD = [scattParam(startIndex + 1, 4),...
        scattParam(startIndex + 2, 2),...
        scattParam(startIndex + 2, 4)];
    for paramIdx = 1:3
        projPot = projPot + 4 * pi^2 * paramA(paramIdx) *...
            besselk(0, 2 * pi * sqrt(rhoSq) * sqrt(paramB(paramIdx))) +...
            2 * pi^2 * paramC(paramIdx) / paramD(paramIdx) *...
            exp(-pi^2 * (rhoSq) / paramD(paramIdx));
    end
end
projPot = bohrRadius * eleCharge * projPot;

end