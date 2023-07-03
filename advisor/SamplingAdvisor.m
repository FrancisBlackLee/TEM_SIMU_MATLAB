function [r] = SamplingAdvisor(keV, Lx, Ly, Nx, Ny, bwl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wavLen = HighEnergyWavLen_X(keV);
r.dx = Lx / Nx;
r.dy = Ly / Ny;
r.dax = wavLen * 1e3 / Lx;
r.day = wavLen * 1e3 / Ly;
r.max = wavLen * 1e3 * bwl / (2 * r.dx);
r.may = wavLen * 1e3 * bwl / (2 * r.dy);

end