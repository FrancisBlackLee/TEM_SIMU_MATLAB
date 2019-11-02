function [MeanDisplace] = ThermoDisplace_0(MassNum, DebyeTemp, Temp)
%ThermoDisplace_0.m calculates the standard deviation of
%thermo-displacement using Einstein model.
%   Input:
%       MassNum -- mass number, i.e. for carbon is 12.01;
%       DebyeTemp -- Debye temperature (in Kelvin);
%       Temp -- simulation temperature;
%   Output:
%       MeanDisplace -- standard deviation of the thermo-displacement (in
%       Angstrom);

MeanDisplace = sqrt(144.38 * Temp ./ (MassNum .* DebyeTemp.^2));

end

