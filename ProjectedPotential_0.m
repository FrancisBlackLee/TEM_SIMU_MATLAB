function Proj_Pot = ProjectedPotential_0(Lx, Ly, Nx, Ny, AtomCoordType)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomZ -- atomic numbers of the input atom series;
%   AtomX, AtomY -- atomic coordinates corresponding to the atomic series;

AtomNum = size(AtomCoordType, 2);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
% delta = 0.1 * dx;
deltaSq = dx * dx + dy * dy;

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom

FileName = mfilename('fullpath');
FileName = strcat(FileName, '.m');
[filepath, name, ext] = fileparts(FileName);
Pot_txt_name = fullfile(filepath, 'Scattering_Factors.txt');
Scatt_Fac = load(Pot_txt_name);

Proj_Pot = zeros(size(X));
for i = 1 : AtomNum
    AtomX = AtomCoordType(1, i);
    AtomY = AtomCoordType(2, i);
    AtomType = AtomCoordType(4, i);
    RHOsq = (X - AtomX).^2 + (Y - AtomY).^2;
    RHOsq(RHOsq < deltaSq) = deltaSq;
    StartIndex = 3 * (AtomType - 1) + 1;
    A = [Scatt_Fac(StartIndex, 1), Scatt_Fac(StartIndex, 3), Scatt_Fac(StartIndex + 1, 1)];
    B = [Scatt_Fac(StartIndex, 2), Scatt_Fac(StartIndex, 4), Scatt_Fac(StartIndex + 1, 2)];
    C = [Scatt_Fac(StartIndex + 1, 3), Scatt_Fac(StartIndex + 2, 1), Scatt_Fac(StartIndex + 2, 3)];
    D = [Scatt_Fac(StartIndex + 1, 4), Scatt_Fac(StartIndex + 2, 2), Scatt_Fac(StartIndex + 2, 4)];
    for j = 1:3
        Proj_Pot = Proj_Pot + 4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(RHOsq) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (RHOsq) / D(j));
    end
end
Proj_Pot = a * e * Proj_Pot;

end

