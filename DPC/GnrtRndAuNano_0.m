% GnrtRndAuNano_0.m -- dynamic diffraction (considering frozen phonon)
clc;
close all;
clear all;
%% Main
% fcc, lattice constant a = 4.0786 Angstrom
LattConst = 4.0786;
IteNum = 21;
vec{1} = [1; 1; 0] / 2;
vec{2} = [0; 1; 1] / 2;
vec{3} = [1; 0; 1] / 2;
growdirec = [-3, -2, -1, 1, 2, 3];

AuNanoSeed = [0; 0; 0];
for IteIdx = 1 : IteNum
    CurPar = AuNanoSeed;
    for AtomIdx = 1 : size(CurPar, 2)
        rng('shuffle');
        grow = rand();
        if grow >= 0.5
            rng('shuffle');
            VecIdx = randperm(6);
            VecIdx = VecIdx(1);
            AuNanoSeed = [AuNanoSeed, CurPar( : , AtomIdx) + sign(growdirec(VecIdx)) * vec{abs(growdirec(VecIdx))}];
        end
    end
    AuNanoSeed = AuNanoSeed';
    AuNanoSeed = uniquetol(AuNanoSeed, 1e-3, 'ByRows', true);
    AuNanoSeed = AuNanoSeed';
end

AuNanoSeed = LattConst * AuNanoSeed;
figure('units','normalized','outerposition',[0 0 1 1]);
scatter3(AuNanoSeed(1, : ), AuNanoSeed(2, : ), AuNanoSeed(3, : ), 'filled');