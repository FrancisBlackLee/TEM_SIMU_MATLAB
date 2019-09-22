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
% SpecimenTF = BandwidthLimit(SpecimenTF, Lx, Ly, Nx, Ny, 0.5);
Params.KeV = 100;
lambda = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
k = 2 * pi / lambda;     %wavenumber
Params.amax = 23;
Params.Cs = 0;
Params.df = 0;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);

%% Scanning Module
Scan_Nx = 128;
Scan_Ny = 128;
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
COM = zeros(2, Scan_Ny * Scan_Nx);
Phase_GradX = zeros(Scan_Nx);
Phase_GradY = zeros(Scan_Ny);
% Add waitbar
W_B = waitbar(0, 'Please wait...');
for i = 1 : Scan_Ny
    yp = ADF_y(i);
    for j = 1 : Scan_Nx
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        Trans_Wave = Probe .* TestTF;
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave)) * dx^2);
        COM(1, (i-1) * Scan_Nx + j) = sum(sum(abs(Trans_Wave_Far.^2) .* Fx))...
                                    / sum(sum(abs(Trans_Wave_Far.^2))); % COM denotes "center of mass"
        COM(2, (i-1) * Scan_Nx + j) = sum(sum(abs(Trans_Wave_Far.^2) .* Fy))...
                                    / sum(sum(abs(Trans_Wave_Far.^2)));
        Phase_GradX(i, j) = 2 * pi * COM(1, (i - 1) * Scan_Nx + j); % Gradients
        Phase_GradY(i, j) = 2 * pi * COM(2, (i - 1) * Scan_Nx + j);
    end
    waitbar(((i - 1) * Scan_Nx + j) / (Scan_Ny * Scan_Nx), W_B, [num2str(((i - 1) * Scan_Nx + j) / (Scan_Ny * Scan_Nx) * 100), '%']);
end
delete(W_B);
%% Retrieve along X and Y respectively
Retri_PotX_L2R = zeros(size(Phase_GradX));
Retri_PotX_L2R( : , 1) = Retri_PotX_L2R( : , 1) + Scan_dx * Phase_GradX( : , 1);
for i = 2 : Scan_Nx
    Retri_PotX_L2R( : , i) = Retri_PotX_L2R( : , i - 1) + Scan_dx * Phase_GradX( : , i);
end

% Retri_PotX_R2L = zeros(size(Phase_GradX));
% Retri_PotX_R2L(:,Scan_N) = Retri_PotX_R2L(:,Scan_N) - Scan_d * Phase_GradX(:,Scan_N);
% for i = Scan_N-1:-1:1
%     Retri_PotX_R2L(:,i) = Retri_PotX_R2L(:,i+1) - Scan_d * Phase_GradX(:,i);
% end
% Retri_PotX_R2L = Retri_PotX_R2L/InterCoeff;

% Retri_PotX = (Retri_PotX_L2R + Retri_PotX_R2L) / 2;
Retri_PotX = Retri_PotX_L2R;
% Retri_PotX = Retri_PotX_R2L;

figure;
subplot(1, 2, 1);
imagesc(ADF_x, ADF_y, Retri_PotX);
colormap('gray');
axis square;
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('DPC Retrieved Potential');
subplot(1, 2, 2);
plot(ADF_x, Retri_PotX(Scan_Ny / 2 + 1, : ), 'b');
axis square;
xlabel('x (in Angstrom)'); ylabel('keV');
title('Potential Profile');

% Retri_PotY = zeros(size(Phase_GradX));
% Retri_PotY(:,1) = Retri_PotY(:,1) + Scan_d * Phase_GradX(:,1);
% for i = 2:Scan_N
%     Retri_PotY(:,i) = Retri_PotY(:,i-1) + Scan_d * Phase_GradX(:,i);
% end
% Retri_PotY = Retri_PotY/InterCoeff;
% figure;
% subplot(1,2,1);
% imagesc(ADF_x,ADF_y,Retri_PotY);
% colormap('gray');
% axis square;
% xlabel('x (in Angstrom)');ylabel('y (in Angstrom)');
% title('DPC Retrieved Potential');
% subplot(1,2,2);
% plot(ADF_x,Retri_PotY(Scan_N/2,:));
% axis square;
% xlabel('x (in Angstrom');ylabel('keV');
% title('Potential Profile');