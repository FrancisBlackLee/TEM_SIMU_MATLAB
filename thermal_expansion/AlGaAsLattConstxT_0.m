function [lattConstAxT] = AlGaAsLattConstxT_0(x, T)
%Calculate lattice constant a for Al_xGa_1-xAs
%   Ref: Adachi, Sadao, ed. Properties of aluminium gallium arsenide. No. 7. IET, 1993.
%   a = a0(1 + ThermalExpanCoeff * DeltaT)

a0 = 5.6533 + 0.0078 * x; % a0 at 300K
ThermalExpanCoeff = (5.97 - 1.76 * x) * 1e-6;
lattConstAxT = a0 .* (1 + ThermalExpanCoeff .* (T - 300));

end

