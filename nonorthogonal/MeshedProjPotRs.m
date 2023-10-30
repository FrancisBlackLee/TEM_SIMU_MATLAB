function projPot = MeshedProjPotRs(xMesh, yMesh, types, occus, axs, ays)
%MeshedProjPotRs calculates projected potential of atoms in real space
%using a pair of predefined x and y meshes.
% projPot = MeshedProjPotRs(xMesh, yMesh, types, occus, axs, ays)
%   xMesh -- mesh of x coordinates;
%   yMesh -- mesh of y coordinates;
%   types -- atomic types;
%   occus -- atomic occupations;
%   axs -- atomic x coordinates;
%   ays -- atomic y coordinates;
%   projPot -- projected potential;

na = length(types);
rSqrThr = 1.0e-8;

a0 = 0.529;
e0 = 14.4;

scattParam = load('Scattering_Factors.txt');

projPot = zeros(size(xMesh));
for ia = 1 : na
    rSqrMesh = (xMesh - axs(ia)).^2 + (yMesh - ays(ia)).^2;
    rSqrMesh(rSqrMesh < rSqrThr) = rSqrThr;
    startIndex = 3 * (types(ia) - 1) + 1;
    A = [scattParam(startIndex, 1), scattParam(startIndex, 3), scattParam(startIndex + 1, 1)];
    B = [scattParam(startIndex, 2), scattParam(startIndex, 4), scattParam(startIndex + 1, 2)];
    C = [scattParam(startIndex + 1, 3), scattParam(startIndex + 2, 1), scattParam(startIndex + 2, 3)];
    D = [scattParam(startIndex + 1, 4), scattParam(startIndex + 2, 2), scattParam(startIndex + 2, 4)];
    for j = 1:3
        projPot = projPot + occus(ia) * 4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(rSqrMesh) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (rSqrMesh) / D(j));
    end
end
projPot = a0 * e0 * projPot;

end