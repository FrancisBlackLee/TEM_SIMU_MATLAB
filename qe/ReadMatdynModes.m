function [qs, bands, eigenVecs] = ReadMatdynModes(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

text = fileread(filename);
% find q = ...
expr = '[^fre]q.*?[\n]';
qInds = regexp(text, expr);
nq = numel(qInds);

% find nBand and nAtom
firstQBlock = text(qInds(1) : qInds(2) - 1);
expr = 'freq.*?\[THz\]';
firstQBlockFreqInds = regexp(firstQBlock, expr);
nBand = numel(firstQBlockFreqInds);
nAtom = nBand / 3;

qs = zeros(nq, 3);
bands = zeros(nq, nBand);
eigenVecs = 1 + 1i * ones(nAtom, 3, nBand, nq);

% parse each q block
for iq = 1 : nq
    if iq ~= nq
        qBlock = text(qInds(iq): qInds(iq + 1) - 1);
    else
        qBlock = text(qInds(iq) : end);
    end

    [qs(iq, :), bands(iq, :), eigenVecs(:, :, :, iq)] = ParseQBlock(qBlock);
end

% nested function:
    function [q, qBands, qEigenVecs] = ParseQBlock(blockStr)

        % match q
        qStr = regexp(qBlock, '[^fre]q.*?[\n]', 'match');
        qStr = regexp(qStr, '[+-]?\d+\.?\d*', 'match');
        q = str2double(qStr{1});

        % match qbands
        bandStrInds = regexp(qBlock, 'freq.*?\[THz\]');
        qBands = zeros(1, nBand);
        qEigenVecs = 1 + 1i * ones(nAtom, 3, nBand);

        for iBand = 1 : nBand
            if iBand ~= nBand
                bandBlock = blockStr(bandStrInds(iBand) : bandStrInds(iBand + 1) - 1);
            else
                bandBlock = blockStr(bandStrInds(iBand) : end);
            end

            [qBands(iBand), qEigenVecs(:, :, iBand)] = ParseBandBlock(bandBlock);
        end
    end


    function [f, vecs] = ParseBandBlock(bandBlock)
        bandStr = regexp(bandBlock, '=.*?\[THz\]', 'match');
        bandStr = regexp(bandStr, '[+-]?\d+\.?\d*', 'match');
        f = str2double(bandStr{1});

        vecStrs = regexp(bandBlock, '[\n].*?\)', 'match');
        vecs = 1 + 1i * ones(nAtom, 3);
        for iAtom = 1 : nAtom
            valStr = regexp(vecStrs{iAtom}, '[+-]?\d+\.?\d*', 'match');
            vals = str2double(valStr);
            vecs(iAtom, :) = vals(1 : 2 : end) + 1i * vals(2 : 2 : end);
        end
    end

end