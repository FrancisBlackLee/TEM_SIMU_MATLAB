%% Graphene
a0 = 1.42;% in Angstrom
a = a0 * sqrt(3);
Nx = 512;
Ny = 512;
Lx = 6 * 3 * a0;
Ly = 10 * a;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
Vz = GrapheneCreate();
SpecimenTF = exp(1i * Vz / 1000);
SpecimenTF = BandwidthLimit(SpecimenTF, Lx, Ly, Nx, Ny, 0.5);
Params.KeV = 200;
lambda = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
k = 2 * pi / lambda;     %wavenumber
Params.amax = 50;
Params.Cs = 0;
Params.df = 0;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);

%% Scanning Module
Scan_Nx = 64;
Scan_Ny = 64;
Scan_Lx = Lx / 2.5;
Scan_Ly = Ly / 2.5;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
ADF_x = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 - Scan_dx;
ADF_y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 - Scan_dy;
[ADF_Xmesh, ADF_Ymesh] = meshgrid(ADF_x, ADF_y);
ADF_FreqX = -1 / (2 * Scan_dx) : 1 / Scan_Lx : 1 / (2 * Scan_dx) - 1 / Scan_Lx;
ADF_FreqY = -1 / (2 * Scan_dy) : 1 / Scan_Ly : 1 / (2 * Scan_dy) - 1 / Scan_Ly;
[ADF_FreqXmesh, ADF_FreqYmesh] = meshgrid(ADF_FreqX, ADF_FreqY);
STEM_Image = zeros(Nx, Ny, Scan_Nx*Scan_Ny);
TestTF = SpecimenTF;
% Show the images of original potential
figure;
subplot(1, 2, 1);
imagesc(x, y, Vz / 1000);
colormap('gray');
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('Original Phase');
subplot(1, 2, 2);
plot(x, Vz(Nx / 2, : ) / 1000);
xlabel('x (in Angstrom'); ylabel('keV');
title('Original phase profile');
for i = 1 : Scan_Ny
    yp = ADF_y(i);
    for j = 1 : Scan_Nx
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        Trans_Wave = Probe .* TestTF;
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave)) * dx * dy);
        STEM_Image( : , : , (i - 1) * Scan_Ny + j) = abs(Trans_Wave_Far .^ 2);
    end
    disp(i / Scan_Ny);
end

%% DPC
% Calculate the center of mass
COM = zeros(1, 2, Scan_Nx * Scan_Ny);
Phase_GradX = zeros(Scan_Ny, Scan_Nx);
Phase_GradY = zeros(Scan_Ny, Scan_Nx);
for i = 1 : Scan_Ny
    for j = 1 : Scan_Nx
        COM(1, 1, (i - 1) * Scan_Ny + j) = sum(sum(STEM_Image( : , : , (i - 1) * Scan_Ny + j) .* Fx))/sum(sum(STEM_Image( : , : , (i - 1) * Scan_Ny + j)));
        COM(1, 2, (i - 1) * Scan_Ny + j) = sum(sum(STEM_Image( : , : , (i - 1) * Scan_Ny + j) .* Fy))/sum(sum(STEM_Image( : , : , (i - 1) * Scan_Ny + j)));
        Phase_GradX(i, j) = 2 * pi * COM(1, 1, (i - 1) * Scan_Ny + j);
        Phase_GradY(i, j) = 2 * pi * COM(1, 2, (i - 1) * Scan_Ny + j);
    end
end
FT_Phase_GradX = fft2(fftshift(Phase_GradX));
temp_X = FT_Phase_GradX ./ (1i * pi * fftshift(ADF_FreqXmesh));
Retri_PhaseX = real(ifftshift(ifft2(temp_X)));
FT_Phase_GradY = fft2(fftshift(Phase_GradY));
temp_Y = FT_Phase_GradY ./ (1i * pi * fftshift(ADF_FreqYmesh));
Retri_PhaseY = real(ifftshift(ifft2(temp_Y)));
Retri_Phase = sqrt(Retri_PhaseX.^2 + Retri_PhaseY.^2);
% Retri_Phase = Retri_PhaseX + Retri_PhaseY;
figure;
subplot(1, 2, 1);
imagesc(ADF_x, ADF_y, Retri_Phase);
colormap('gray');
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('DPC Retrieved Phase');
subplot(1, 2, 2);
plot(ADF_x, unwrap(Retri_Phase(Scan_Ny / 2, : )));
xlabel('x (in Angstrom'); ylabel('keV');
title('Phase Profile');