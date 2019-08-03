clc
close all
clear all

%% basic information
L = 30;     %side length 
M = 1024;      %number of samples
dx = L / M;   %scr sample interval
x = -L / 2 : dx : L / 2 -dx;  %scr coords
y = x;
Params.KeV = 200;
lambda = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
k = 2 * pi / lambda;     %wavenumber
InterCoeff = InteractionCoefficient(Params.KeV);
Params.amax = 23;
Params.Cs = 0;
Params.df = 0;

[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / L : 1 / (2 * dx) - 1 / L;fy = fx;
[Fx, Fy] = meshgrid(fx, fy);

% Draw probe
Probe_test = ProbeCreate(Params, 0, 0, L, L, M, M);
Probe_test_Intensity = abs(Probe_test.^2);
figure;
subplot(1, 2, 1);
imagesc(x, y, Probe_test_Intensity);
colormap('gray');
title('Probe');
subplot(1, 2, 2);
mesh(x, y, Probe_test_Intensity);

%% Table Example
R = sqrt((X - 0 * ones(size(X))).^2 + (Y - 0 * ones(size(Y))).^2);
Table_Potential = -1000 * R + 2000 * ones(size(R));
Table_Potential(R <= 0.5) = 1500;
Table_Potential(R >= 2) = 0;
TableTF = exp(InterCoeff * 1i * Table_Potential / 1000);
% TableTF = BandwidthLimit(TableTF,L,L,M,M,0.5);

%% Single atom ADF-STEM
Proj_Pot = ProjectedPotential(L, L, M, M, 6, -8, 0) + ...
           ProjectedPotential(L, L, M, M, 14, -4, 0) + ...
           ProjectedPotential(L, L, M, M, 29, 0, 0) + ...
           ProjectedPotential(L, L, M, M, 79, 4, 0) + ...
           ProjectedPotential(L, L, M, M, 92, 8, 0);
SpecimenTF = exp(InterCoeff * 1i * Proj_Pot / 1000);   %Proj_ is the projected atomic potential
% SpecimenTF = BandwidthLimit(SpecimenTF,L,L,M,M,0.5);

%% Scanning module
% Scanning parameters:
Scan_N = 256;
Scan_L = L / 1.5;
Scan_d = Scan_L / Scan_N;
ADF_x = -Scan_L / 2 : Scan_d : Scan_L / 2 - Scan_d;
ADF_y = ADF_x;
[ADF_Xmesh, ADF_Ymesh] = meshgrid(ADF_x, ADF_y);
ADF_FreqX = -1 / (2 * Scan_d) : 1 / Scan_L : 1 / (2 * Scan_d) - 1 / Scan_L;
ADF_FreqY = ADF_FreqX;
[ADF_FreqXmesh, ADF_FreqYmesh] = meshgrid(ADF_FreqX, ADF_FreqY);

COM = zeros(2, Scan_N * Scan_N);
Phase_GradX = zeros(Scan_N);
Phase_GradY = zeros(Scan_N);
TestTF = SpecimenTF; % Choose transfer function
% Show the images of original potential
Test_Potential = Proj_Pot; % Choose potential to image and plot
figure;
subplot(1, 2, 1);
imagesc(x, y, Test_Potential / 1000);
colormap('gray');
axis square;
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('Original Potential');
subplot(1, 2, 2);
plot(x, Test_Potential(M/2, : ) / 1000, 'b');
axis square;
xlabel('x (in Angstrom)'); ylabel('keV');
title('Original Potential profile');

% Add waitbar
W_B = waitbar(0, 'Please wait...');
for i = 1 : Scan_N
    yp = ADF_y(i);
    for j = 1 : Scan_N
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, L, L, M, M);
        Trans_Wave = Probe .* TestTF;
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave)) * dx^2);
        COM(1, (i-1) * Scan_N + j) = sum(sum(abs(Trans_Wave_Far.^2) .* Fx))...
                                    / sum(sum(abs(Trans_Wave_Far.^2))); % COM denotes "center of mass"
        COM(2, (i-1) * Scan_N + j) = sum(sum(abs(Trans_Wave_Far.^2) .* Fy))...
                                    / sum(sum(abs(Trans_Wave_Far.^2)));
        Phase_GradX(i, j) = 2 * pi * COM(1, (i - 1) * Scan_N + j); % Gradients
        Phase_GradY(i, j) = 2 * pi * COM(2, (i - 1) * Scan_N + j);
    end
    waitbar(((i - 1) * Scan_N + j) / (Scan_N * Scan_N), W_B, [num2str(((i - 1) * Scan_N + j) / (Scan_N * Scan_N) * 100), '%']);
end
delete(W_B);

%% DPC Electric Field Retrieval
% The original electric field
[EleField_x, EleField_y] = gradient(Test_Potential, dx, dx);
Ori_EleField = wasted_temp();
% Retrieved electric field
ElectricField = sqrt(Phase_GradX.^2 + Phase_GradY.^2);
% figure;
% imagesc(ADF_x, ADF_y, ElectricField);
% colormap('gray');
% axis square;
% xlabel('x (in Angstrom)');ylabel('y (in Angstrom)');
% title('Electric Field');
figure;
subplot(1, 2, 1);
contour(x, y, Test_Potential);
hold on;
quiver(x, y, -EleField_x, -EleField_y);
hold off;
axis square; 
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('Original Electric Field');
subplot(1, 2, 2);
contour(ADF_x, ADF_y, ElectricField);
hold on;
quiver(ADF_x, ADF_y, -Phase_GradX, -Phase_GradY);
hold off;
axis square;
xlabel('x (in Angstrom)'); ylabel('y (in Angstrom)');
title('Retrived Electric Field');
figure;
plot(ADF_x, Ori_EleField(Scan_N/2, :), 'r', ADF_x, ElectricField(Scan_N/2, :), 'b');
legend('Original', 'Retrieved');
xlabel('x (Angstrom)');
ylabel('kV / Ang.');
title('Electric field');

%% Retrieve along X and Y respectively
Retri_PotX_L2R = zeros(size(Phase_GradX));
Retri_PotX_L2R( : , 1) = Retri_PotX_L2R( : , 1) + Scan_d * Phase_GradX( : , 1);
for i = 2 : Scan_N
    Retri_PotX_L2R( : , i) = Retri_PotX_L2R( : , i - 1) + Scan_d * Phase_GradX( : , i);
end
Retri_PotX_L2R = Retri_PotX_L2R / InterCoeff;

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
% subplot(1,2,1);
% imagesc(ADF_x,ADF_y,Retri_PotX);
% colormap('gray');
% axis square;
% xlabel('x (in Angstrom)');ylabel('y (in Angstrom)');
% title('DPC Retrieved Potential');
% subplot(1,2,2);
plot(ADF_x, Retri_PotX(Scan_N / 2, : ), 'b');
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