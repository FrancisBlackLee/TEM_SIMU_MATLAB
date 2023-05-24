function [pwscf] = ExpandPwscfAtomMass(pwscf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pwscf.mass = zeros(1, pwscf.system.nat);
for iType = 1 : pwscf.system.ntyp
    pwscf.mass(pwscf.atomic_positions.types == pwscf.atomic_species.types(iType)) =...
        pwscf.atomic_species.masses(iType);
end

end