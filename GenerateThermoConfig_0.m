function [] = GenerateThermoConfig_0(SavePath, StdLatt, MassNum, DebyeTemp, SimuTemp, ConfigNum)
%GenerateThermoConfig_0.m generates a given number of lattice
%configurations under input temperature using Einstein model.
%   SavePath -- save all the generated configuration under a given path, to
%       avoid memory overflow;
%   StdLatt -- standard lattice matrix without displacement, matrix format
%       is [T1, ..., TN; P1, ..., PN; X1, ..., XN; Y1, ..., YN; Z1, ...,
%       ZN], where T denotes atomic type, represented by the atomic
%       numbers, P denotes the elemental proportion, X, Y and Z denote the
%       orthogonal coordinates;
%   MassNum -- mass number array corresponding to the lattice matrix;
%   DebyeTemp -- Debye temperature array corresponding to the lattice
%       matrix;
%   SimuTemp -- simulation temperature;
%   ConfigNum -- number of configurations to be generated;

ThermoDisp = ThermoDisplace_0(MassNum, DebyeTemp, SimuTemp);

for ConfigIdx = 1 : ConfigNum
    % normally distributed random radial displacement:
    rng('shuffle');
    RndRadius = normrnd(zeros(size(ThermoDisp)), ThermoDisp);
    % uniformly distributed random azimuthal angles:
    rng('shuffle');
    RndAzi = pi * rand(size(ThermoDisp));
    % uniformly distributed random polar angles:
    rng('shuffle');
    RndPol = 2 * pi * rand(size(ThermoDisp));

    Xdisp = RndRadius .* sin(RndAzi) .* cos(RndPol);
    Ydisp = RndRadius .* sin(RndAzi) .* sin(RndPol);
    Zdisp = RndRadius .* cos(RndAzi);
    
    DispLatt = StdLatt;
    DispLatt(3, : ) = DispLatt(3, : ) + Xdisp;
    DispLatt(4, : ) = DispLatt(4, : ) + Ydisp;
    DispLatt(5, : ) = DispLatt(5, : ) + Zdisp;
    
    DispLatt = DispLatt';
    filename = strcat(SavePath, 'Config_', num2str(ConfigIdx), '.txt');
    save(filename, 'DispLatt', '-ascii', '-double', '-tabs');
end

end

