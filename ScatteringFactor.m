function [sf] = ScatteringFactor(AtomType, q)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FileName = mfilename('fullpath');
FileName = strcat(FileName, '.m');
[filepath, name, ext] = fileparts(FileName);
Pot_txt_name = fullfile(filepath, 'Scattering_Factors.txt');
Scatt_Fac = load(Pot_txt_name);

StartIndex = 3 * (AtomType - 1) + 1;
A = [Scatt_Fac(StartIndex, 1), Scatt_Fac(StartIndex, 3), Scatt_Fac(StartIndex + 1, 1)];
B = [Scatt_Fac(StartIndex, 2), Scatt_Fac(StartIndex, 4), Scatt_Fac(StartIndex + 1, 2)];
C = [Scatt_Fac(StartIndex + 1, 3), Scatt_Fac(StartIndex + 2, 1), Scatt_Fac(StartIndex + 2, 3)];
D = [Scatt_Fac(StartIndex + 1, 4), Scatt_Fac(StartIndex + 2, 2), Scatt_Fac(StartIndex + 2, 4)];

sf = 0;
for i = 1 : 3
    sf = sf + A(i)/(q^2 + B(i)) + C(i)*exp(-D(i)*q^2);
end
end

