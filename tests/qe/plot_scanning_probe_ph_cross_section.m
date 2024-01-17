% plot_scanning_probe_ph_cross_section.m
clc;
clear;
close all;
%% main:
nProbe = 5;
dmeV = 0.1;
meV = 0 : dmeV : 80;
nmeV = length(meV);
phSpectra = zeros(1, nmeV);

colors = ["k", "r", "g", "y", "b"];

figure;
hold on;
for iProbe = 1 : nProbe
    filename = ['eloss_data_probe_', num2str(iProbe), '.mat'];
    load(fullfile('tests/qe', filename));

    for imeV = 1 : nmeV
        phSpectra(imeV) = sum(phCs(eLosses > meV(imeV) - 0.5 * dmeV & ...
            eLosses < meV(imeV) + 0.5 * dmeV), 'all');
    end
    
    filtPhSpectra = gaussfilt(meV, phSpectra, 1);

%     plot(meV, filtPhSpectra, 'LineWidth', 2);
    plot3(meV, iProbe * ones(size(meV)), filtPhSpectra, 'LineWidth', 2, 'Color', colors(iProbe));
    fill3(meV, iProbe * ones(size(meV)), filtPhSpectra, colors(iProbe));
%     lgdStr{iProbe} = ['probe ', num2str(iProbe)];
end
hold off;
xlabel('energy loss (meV)');
ylabel('probe');
zlabel('cross section (a.u.)');
% legend(lgdStr);
axis square;