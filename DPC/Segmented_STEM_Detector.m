function [Detector_array, COM] = Segmented_STEM_Detector(DetectorAngles, Seg_N, Lx, Ly, Nx, Ny, lambda, ShowDetectorYN)
%UNTITLED2 Summary of this function goes here
%   Input:
%       DetectorParams --detector parameters, usually about its shape
%       Seg_N --number of segments.
%       Lx, Ly --side length
%       Nx, Ny --number of sampling
%       lambda --wavelength
%       ShowDetectorYN --whether to show the images of the detector
%   Output:
%       Detector_array --array of detectors

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
fx = linspace(-1/(2*dx),1/(2*dx),Nx);
fy = linspace(-1/(2*dy),1/(2*dy),Ny);
[Fx,Fy] = meshgrid(fx,fy);
FreqSquare = Fx.^2+Fy.^2;
switch Seg_N
    case 2
        D_1 = ones(size(Fx));
        D_2 = ones(size(Fx));
        D_1((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx>0)) = 0;
        D_2((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx<0)) = 0;
        Detector_array = cat(3,D_1,D_2);
        k1 = tan(DetectorAngles(1)*0.001)/lambda;
        k2 = tan(DetectorAngles(2)*0.001)/lambda;
        k_COM = 4/(3*pi)*(k2^3-k1^3)/(k2^2-k1^2);
        COM = [-k_COM, k_COM];
    case 4
        D_1 = ones(size(Fx));
        D_2 = ones(size(Fx));
        D_3 = ones(size(Fx));
        D_4 = ones(size(Fx));
        D_1((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx<0)|(Fy>0)) = 0;
        D_2((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx<0)|(Fy<0)) = 0;
        D_3((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx>0)|(Fy<0)) = 0;
        D_4((FreqSquare<(tan(DetectorAngles(1)*0.001)/lambda)^2)|(FreqSquare>(tan(DetectorAngles(2)*0.001)/lambda)^2)|(Fx>0)|(Fy>0)) = 0;
        D_1 = imrotate(D_1,-45,'bilinear','crop');
        D_2 = imrotate(D_2,-45,'bilinear','crop');
        D_3 = imrotate(D_3,-45,'bilinear','crop');
        D_4 = imrotate(D_4,-45,'bilinear','crop');
        Detector_array = cat(3,D_1,D_2,D_3,D_4);
        k1 = tan(DetectorAngles(1)*0.001)/lambda;
        k2 = tan(DetectorAngles(2)*0.001)/lambda;
        k_COM = 4*sqrt(2)/(3*pi)*(k2^3-k1^3)/(k2^2-k1^2);
        COM = [k_COM, 1i*k_COM, -k_COM, -1i*k_COM];
    otherwise
        disp('This type of detector is not included');
end
if ShowDetectorYN == 1
    for i = 1:Seg_N/2
        subplot(Seg_N/2,2,2*i-1);
        imagesc(Detector_array(:,:,2*i-1));
        colormap('gray');
        title(['Segment',num2str(2*i-1)]);
        subplot(Seg_N/2,2,2*i);
        imagesc(Detector_array(:,:,2*i));
        colormap('gray');
        title(['Segment',num2str(2*i)]);
    end
end

end

