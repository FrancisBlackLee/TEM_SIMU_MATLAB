function [InterCoeff] = InteractionCoefficient(KeV)
%InteractionCofficient.m
%   KeV --electron beam energy (in kilo-electron-volt).

LightSpeed = 2.998e8; % in meter
EleCharge = 1.602e-19; % in coulomb
V = 1000 * KeV; % in volt
m0 = 9.109e-31; % in kilogram
lambda = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);
InterCoeff = 2 * pi/(lambda * KeV) * (m0 * LightSpeed^2 + EleCharge * V)/(2 * m0 * LightSpeed^2 + EleCharge * V);

end

