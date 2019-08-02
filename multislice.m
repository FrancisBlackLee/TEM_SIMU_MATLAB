function [ExitWave] = multislice(IncidentWave, WaveLength, Lx, Ly, TransFuncs, LayerDist, StackNum)
%MULTISLICE.M performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy.
%   IncidentWave;
%   Lx, Ly -- sampling parameters;
%   TransFuncs -- transmission functions;
%   LayerDist -- distance from one layer to the next layer;
%   StackNum -- number of stackings.

Wave = IncidentWave;
LayerNum = length(LayerDist);
for i = 1: StackNum
    for j = 1: LayerNum
        Wave = Wave .* TransFuncs(:, :, j);
        Wave = propTF_1(Wave, Lx, Ly, WaveLength, LayerDist(j));
    end
end
ExitWave = Wave;

end

