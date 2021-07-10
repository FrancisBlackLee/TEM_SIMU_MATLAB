function [expanCoeff] = SiliconThermalExpan_0(T)
%SiliconThermalExpan_0.m calculates the thermal expansion coefficient.
%   Ref: Yim, W.M., and R.J. Paff, J. Appl. Phys. 45, (1974) 1456-1457.

aT0 = 5.4304;
aT = aT0 + 1.8138e-5 * T + 1.542e-9 * T.^2;
expanCoeff = aT / aT0;

end