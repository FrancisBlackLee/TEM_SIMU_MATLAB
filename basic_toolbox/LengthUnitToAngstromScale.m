function scaleCoeff = LengthUnitToAngstromScale(lengthUnit)
%LengthUnitToAngstromScale() generates the scaling coefficient for
%transforming input length unit to angstrom.
% Input:
%   lengthUnit -- valid input: mm, um, nm, pm
% Output:
%   scaleCoeff -- scaling coefficient, 0 for invalid input

switch lengthUnit
    case 'mm'
        scaleCoeff = 1e7;
    case 'um'
        scaleCoeff = 1e4;
    case 'nm'
        scaleCoeff = 10;
    case 'pm'
        scaleCoeff = 1e-2;
    otherwise
        scaleCoeff = 0;
end

end

