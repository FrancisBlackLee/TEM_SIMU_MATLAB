function [seAmp] = SecondaryAmplitude(keV, atomType, atomX, atomY, Lx, Ly, Nx, Ny)
%SecondaryAmplitude.m calculates the secondary electron amplitude for a
%single atom. Ref: Wu L, Egerton R F, Zhu Y. Image simulation for atomic 
%resolution secondary electron image[J]. Ultramicroscopy, 2012, 123: 66-73.
%   keV -- energy of incident beam (in keV);
%   atomType -- atom type represented by it atomic number;
%   atomX -- x coordinate of the atom;
%   atomY -- y coordinate of the atom;
%   Lx, Ly -- side lengths;
%   Nx, Ny -- number of pixels in x and y dimension;

energyLoss = 6.75 * atomType / 1000;
energyLossAngle = (keV + 511) / (keV + 1022) * (energyLoss / keV);
energyLossAngleSqr = energyLossAngle^2;
rolloffAngleSqr = 2 * energyLossAngle;

wavLen = HighEnergyWavLen_X(keV);
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[fxMesh, fyMesh] = meshgrid(fx, fy);
angleSqrMesh = wavLen^2 * (fxMesh.^2 + fyMesh.^2);

diffCrossSect = (energyLossAngleSqr ./ (energyLossAngleSqr + angleSqrMesh)) .*...
    (1 ./ (angleSqrMesh + rolloffAngleSqr));

seAmp = sqrt(diffCrossSect) .* exp(-1i * 2 * pi * (fxMesh * atomX + fyMesh * atomY));
seAmp = 4 * pi^2 * abs(ifftshift(ifft2(fftshift(seAmp))));

end

