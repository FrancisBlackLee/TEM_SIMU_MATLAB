function v = IsoAtomMeanPotential(vol, atypes, aoccus)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

planckConst = 6.62607015e-34; % m^2 * kg / s
eMass = 9.1093837015e-31; % kg
e0 = 1.60217663e-19; % coulomb
v = 0;
for ia = 1 : length(atypes)
    v = v + aoccus(ia) * ScatteringFactor(atypes(ia), 0);
end
v = planckConst^2 / (2 * pi * eMass * e0 * vol * 1.0e-30) * v * 1e-10;

end