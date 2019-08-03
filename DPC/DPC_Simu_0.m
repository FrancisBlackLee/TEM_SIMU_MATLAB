clc
close all
clear all

%% basic information
L = 10;     %side length 
M = 256;      %number of samples
dx = L / M;   %scr sample interval
x = -L / 2 : dx : L / 2 -dx;  %scr coords
y = x;
Params.KeV = 200;
lambda = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
k = 2 * pi / lambda;     %wavenumber
Params.amax = 10.37;
Params.Cs = 1.3;
Params.df = 700;

[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / L : 1 / (2 * dx) - 1 / L;
fy = fx;
[Fx, Fy] = meshgrid(fx,fy);

% Parameters of the detector:
Seg_N = 4; % Note that Seg_N is the number of detector segments, so are the following Seg_N's.
DetectorAngles = [0, 100]; % Semiangle of the detector (in mrad).
[Detector_array, COM] = Segmented_STEM_Detector(DetectorAngles, Seg_N, L, L, M, M, lambda, 1);
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
R = sqrt(X.^2 + Y.^2);
Table_Potential = -250 * R + 500 * ones(size(R));
% Table_Potential(R<=1) = 500;
Table_Potential(R >= 2) = 0;
figure;
plot(x, Table_Potential(M / 2 + 1, : ));
title('Original phase profile');
TableTF = exp(1i * Table_Potential / 1000);
% TableTF = BandwidthLimit(TableTF,L,L,M,M,0.67);

%% Scanning module
% Scanning parameters:
Scan_N = 32;
Scan_L = L / 2;
Scan_d = Scan_L / Scan_N;
ADF_x = -Scan_L / 2 : Scan_d : Scan_L / 2 - Scan_d;
ADF_y = ADF_x;
ADF_FreqX = -1 / (2 * Scan_d) : 1 / Scan_L : 1 / (2 * Scan_d) - 1 / Scan_L;
ADF_FreqY = ADF_FreqX;
[ADF_FreqXmesh, ADF_FreqYmesh] = meshgrid(ADF_FreqX, ADF_FreqY);
STEM_Image = zeros(Scan_N, Scan_N, Seg_N);
TestTF = TableTF;
for i = 1 : Scan_N
    yp = ADF_y(i);
    for j = 1 : Scan_N
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, L, L, M, M);
        Trans_Wave = Probe .* TestTF;
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave)) * dx^2);
        for Det_Seg_index = 1 : Seg_N
            Det_Signal( : , : , Det_Seg_index) = abs(Trans_Wave_Far.^2) .* Detector_array( : , : , Det_Seg_index);
            STEM_Image(i, j, Det_Seg_index) = sum(sum(Det_Signal( : , : , Det_Seg_index)));
        end
    end
    disp(i / Scan_N);
end
% Show the STEM images produced by each detector
figure;
for i = 1 : Seg_N/2
    subplot(Seg_N/2, 2, 2 * i - 1);
    imagesc(ADF_x, ADF_y, STEM_Image( : , : , 2 * i - 1));
    colormap('gray');
    title(['image produced by', 'Segment', num2str(2 * i - 1)]);
    subplot(Seg_N / 2, 2, 2 * i);
    imagesc(ADF_x, ADF_y, STEM_Image( : , : , 2 * i));
    colormap('gray');
    title(['image produced by', 'Segment', num2str(2 * i)]);
end
%% DPC

Intensity_R = zeros(Scan_N);
for i = 1 : Seg_N
    Intensity_R = Intensity_R + COM(i) * STEM_Image( : , : , i);
end
Probe_R = ProbeCreate(Params, 0, 0, Scan_L, Scan_L, Scan_N, Scan_N);
Probe_R_Sq = abs(Probe_R.^2);
Probe_R_Sq(Probe_R_Sq < 1e-10) = 1e-10;
FreqCoords_R = ADF_FreqXmesh + 1i * ADF_FreqYmesh;
FreqCoords_R(abs(FreqCoords_R.^2) < 1e-10) = 1e-10;
Reci_nume = fft2(fftshift(Intensity_R));
Reci_deno = 1i * fftshift(FreqCoords_R) .* fft2(fftshift(Probe_R_Sq));
Retri_Phase = ifftshift(ifft2(Reci_nume ./ Reci_deno));
figure;
subplot(1, 2, 1);
imagesc(ADF_x, ADF_y, abs(Retri_Phase.^2));
colormap('gray');
title('DPC retrieved phase');
subplot(1, 2, 2);
hold on;
plot(ADF_x, Retri_Phase(Scan_N / 2 + 1, : ));
title('Profile');