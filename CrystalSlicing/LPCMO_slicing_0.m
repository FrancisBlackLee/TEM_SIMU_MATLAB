% LPCMO_slicing_0.m
clc;
clear all;
close all;
%% Find the new crystal constants:
W_B = waitbar(0, 'Preparing the specimen...');

FileDir = 'D:\Francis. B. Lee\cooperation\Group Cooperation\Wei_Cao\';
CrysMat = load(strcat(FileDir, 'LMO_101_Coordp.txt'));
CrysMat = [CrysMat( : , 1), ones(size(CrysMat( : , 1))), CrysMat( : , 6), CrysMat( : , 7), CrysMat( : , 8)]';

% Fill Pr and Ca at the positions of La, and thus adapting its elemental
% fraction:
LMO_AtomNum = size(CrysMat, 2);

for i = 1 : LMO_AtomNum
    if CrysMat(1, i) == 57
        CrysMat(2, i) = 0.425;
        CrysMat = [CrysMat, [59; 0.2; CrysMat(3 : 5, i)]];
        CrysMat = [CrysMat, [20; 0.375; CrysMat(3 : 5, i)]];
    end
end

[Z, Order] = sort(CrysMat(5, : ));
TestMat = CrysMat( : , Order)';
ConstA = 5.4430;
ConstB = 7.6855;
ConstC = 5.4490;
ConstAp = 2 * ConstA * ConstC / sqrt(ConstA^2 + ConstC^2);
ConstBp = ConstB;
ConstCp = sqrt(ConstA^2 + ConstC^2);
ConstACp = (ConstAp + ConstCp) / 2;
CrysMat(3, : ) = CrysMat(3, : ) + ConstAp / 2;
CrysMat(5, : ) = CrysMat(5, : ) + ConstCp;
FracCrysMat = CrysMat;
FracCrysMat(3, : ) = FracCrysMat(3, : ) / ConstAp;
FracCrysMat(4, : ) = FracCrysMat(4, : ) / ConstBp;
FracCrysMat(5, : ) = FracCrysMat(5, : ) / ConstCp;

DistError = 1;
UnitIdx = find((FracCrysMat(3, : ) >= 0) & (FracCrysMat(3, : ) <= 1)...
             & (FracCrysMat(4, : ) >= 0) & (FracCrysMat(4, : ) <= 1)...
             & (FracCrysMat(5, : ) >= 0) & (FracCrysMat(5, : ) <= 1));
FracCrysMat = FracCrysMat( : , UnitIdx);
LattConst = [ConstACp, ConstBp, ConstACp];
CellNum = [4, 4, 1];

% Expand the lattice
ExpanCrysMat = SquareLattExpanX(FracCrysMat, LattConst, CellNum);
ExpanAtomNum = size(ExpanCrysMat, 2);
PlotColor = ones(1, size(ExpanCrysMat, 2));
for i = 1 : ExpanAtomNum
    if ExpanCrysMat(1, i) == 25
        PlotColor(i) = 2;
    elseif ExpanCrysMat(1, i) == 8
        PlotColor(i) = 3;
    end
end

[slice, SliceDist] = CrystalSlicing_X(ExpanCrysMat, DistError, LattConst(3), 0, PlotColor);
%% Sampling setting:
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% STEM settings:
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
Params.amax = 22;
Params.Cs = 0;
Params.df = 0;

detector = Fx.^2 + Fy.^2;
detector_cri = detector;
HighAngle = 200 * 0.001;
LowAngle = 40 * 0.001;
detector((detector_cri > (sin(LowAngle) / WaveLength)^2) & (detector_cri < (sin(HighAngle) / WaveLength)^2)) = 1;
detector((detector_cri < (sin(LowAngle) / WaveLength)^2) | (detector_cri > (sin(HighAngle) / WaveLength)^2)) = 0;
%% Generate projected potential:
for i = 1 : length(SliceDist)
    PotFileName = strcat(FileDir, 'SliceProjPot_', num2str(i), '.txt');
    
    ProjPot = ProjectedPotential_1(Lx, Ly, Nx, Ny, slice{i});
    save(PotFileName, 'ProjPot', '-ascii', '-double', '-tabs');
%     % show the projected potential of each slice;
    figure;
    imagesc(x, y, ProjPot);
    colormap('gray'); axis square; title(['z = ', num2str(SliceDist(i))]);
    
    SliceTF = exp(1i * InterCoeff * ProjPot / 1000);
    SliceTF = BandwidthLimit(SliceTF, Lx, Ly, Nx, Ny, 2/3);
    % Put these transmission functions into a 3d array
    TransFuncs(:, :, i) = SliceTF;
%     % Check the transmission functions:
%     figure; subplot(1, 2, 1);
%     imagesc(x, y, real(SliceTF)); colormap('gray'); axis square;
%     title('real'); subplot(1, 2, 2); imagesc(x, y, imag(SliceTF));
%     colormap('gray'); axis square; title('imag');
end

waitbar(0, W_B, 'Specimen preparation completed, start scanning...');
%% Scanning module
% Scanning parameters:
Scan_Nx = 128; % scanning sampling number, adaptive
Scan_Ny = 128;
Scan_Lx = 25.6; % scanning side length, adaptive
Scan_Ly = 25.6;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
ADF_x = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 - Scan_dx;
ADF_y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 - Scan_dy;
STEM_IMAGE = zeros(Scan_Ny, Scan_Nx);
figure;
StackNum = 80; % determines the thickness of the specimen
TotalNum = Scan_Ny * Scan_Nx;
for i=1:Scan_Ny
    yp = ADF_y(i);
    for j=1:Scan_Nx
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        Trans_Wave = multislice(Probe, WaveLength, Lx, Ly, TransFuncs, SliceDist, StackNum);
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave))*dx^2);
        DetectInten = abs(Trans_Wave_Far.^2).*detector;
        STEM_IMAGE(i,j) = sum(sum(DetectInten));
        CurrentNum = (i - 1) * Scan_Nx + j;
        imagesc(ADF_x, -ADF_y, STEM_IMAGE);
        map = colormap(gray);
        axis square;
        title('Example');
        drawnow;
        % Save gif to local directory
%         F = getframe(gcf); 
%         I = frame2im(F); 
%         [I, map] = rgb2ind(I, 256); 
%         if CurrentNum == 1 
%             imwrite(I,map,'D:\Francis. B. Lee\Practice\Conventional Multislice in MATLAB\Specimen_Thickness\Secret\secret.gif','gif','Loopcount',inf,'DelayTime',0.02); 
%         else
%             imwrite(I,map,'D:\Francis. B. Lee\Practice\Conventional Multislice in MATLAB\Specimen_Thickness\Secret\secret.gif','gif','WriteMode','append','DelayTime',0.02); 
%         end
        waitbar(roundn(CurrentNum / TotalNum, -3), W_B, [num2str(roundn((CurrentNum / TotalNum), -3) * 100), '%']);
    end
end
delete(W_B);